"""Build a GitHub dependency submission snapshot from ``uv.lock``.

GitHub's static ``uv.lock`` parser reports every package as ``runtime``, so dev-only packages
(pytest, python-semantic-release -> gitpython, ...) would fail ``actions/dependency-review-action``
on advisories that never reach users. User submissions take priority over static parsing for the
same manifest, so we submit the ``uv.lock`` manifest ourselves with accurate scopes:

* ``runtime``     - the package is reachable from the project's ``dependencies`` or
  ``optional-dependencies`` (installed for users of the package), directly or transitively,
* ``development`` - it is only reachable from dependency groups (``dev``, ``release``, ...).

Usage (the snapshot is written to stdout)::

    python .github/tools/uv_dependency_snapshot.py --sha SHA --ref REF --job-id ID | \\
        gh api -X POST repos/OWNER/REPO/dependency-graph/snapshots --input -

Docs: https://docs.github.com/en/rest/dependency-graph/dependency-submission
"""

import argparse
import json
import os
import re
import sys
import tomllib
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path

MANIFEST = "uv.lock"
CORRELATOR = "uv-lock-scopes"
DETECTOR = {
    "name": "uv-lock-scopes",
    "version": "1.0.0",
    "url": "{server}/{repo}/blob/master/.github/tools/uv_dependency_snapshot.py".format(
        server=os.environ.get("GITHUB_SERVER_URL", "https://github.com"),
        repo=os.environ.get("GITHUB_REPOSITORY", "NaturalAntibody"),
    ),
}

Key = tuple[str, str]  # (normalized name, version)


def normalize(name: str) -> str:
    """PEP 503 name normalization."""
    return re.sub(r"[-_.]+", "-", name).lower()


def purl(name: str, version: str) -> str:
    return f"pkg:pypi/{name}@{version}"


def is_project(package: dict) -> bool:
    """The project itself (and workspace members) are editable/virtual local sources."""
    source = package.get("source", {})
    return "editable" in source or "virtual" in source


class Lock:
    def __init__(self, lock: dict) -> None:
        self.packages: dict[Key, dict] = {}
        self.by_name: dict[str, list[Key]] = defaultdict(list)
        for package in lock["package"]:
            if "version" not in package:  # dynamic-version local packages
                continue
            key = (normalize(package["name"]), package["version"])
            self.packages[key] = package
            self.by_name[key[0]].append(key)

    def resolve(self, dep: dict) -> list[Key]:
        """Lock entries a dependency edge points to (all versions when the edge does not pin one)."""
        name = normalize(dep["name"])
        if "version" in dep:
            return [(name, dep["version"])]
        return list(self.by_name.get(name, []))

    def edges(self, key: Key, extras: frozenset[str]) -> list[dict]:
        package = self.packages[key]
        deps = list(package.get("dependencies", []))
        for extra in extras:
            deps += package.get("optional-dependencies", {}).get(extra, [])
        return deps


def runtime_closure(lock: Lock, roots: list[dict]) -> set[Key]:
    """Every package reachable from the given dependency edges, following requested extras."""
    seen: set[tuple[Key, frozenset[str]]] = set()
    queue = deque(roots)
    while queue:
        dep = queue.popleft()
        extras = frozenset(dep.get("extra", []))
        for key in lock.resolve(dep):
            if key not in lock.packages or (key, extras) in seen:
                continue
            seen.add((key, extras))
            queue.extend(lock.edges(key, extras))
    return {key for key, _ in seen}


def build_manifest(raw_lock: dict) -> dict:
    lock = Lock(raw_lock)
    projects = [p for p in lock.packages.values() if is_project(p)]

    runtime_roots: list[dict] = []
    direct: set[str] = set()
    for project in projects:
        runtime_roots += project.get("dependencies", [])
        for deps in project.get("optional-dependencies", {}).values():
            runtime_roots += deps
        direct.update(normalize(d["name"]) for d in runtime_roots)
        for deps in project.get("dev-dependencies", {}).values():
            direct.update(normalize(d["name"]) for d in deps)
    runtime = runtime_closure(lock, runtime_roots)

    resolved = {}
    for key, package in lock.packages.items():
        if is_project(package):
            continue
        name, version = key
        children = sorted(
            purl(*child)
            for dep in lock.edges(key, frozenset(package.get("optional-dependencies", {})))
            for child in lock.resolve(dep)
            if child in lock.packages and not is_project(lock.packages[child])
        )
        package_url = purl(name, version)
        resolved[package_url] = {
            "package_url": package_url,
            "relationship": "direct" if name in direct else "indirect",
            "scope": "runtime" if key in runtime else "development",
            "dependencies": children,
        }
    return {"name": MANIFEST, "file": {"source_location": MANIFEST}, "resolved": resolved}


def build_snapshot(lock: dict, sha: str, ref: str, job_id: str) -> dict:
    return {
        "version": 0,
        "sha": sha,
        "ref": ref,
        "job": {"correlator": CORRELATOR, "id": job_id},
        "detector": DETECTOR,
        "scanned": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "manifests": {MANIFEST: build_manifest(lock)},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sha", required=True, help="commit SHA the snapshot belongs to")
    parser.add_argument("--ref", required=True, help="full git ref, e.g. refs/heads/master")
    parser.add_argument("--job-id", required=True, help="CI run id")
    parser.add_argument("--lock", type=Path, default=Path(MANIFEST))
    args = parser.parse_args()

    with args.lock.open("rb") as file:
        lock = tomllib.load(file)

    json.dump(build_snapshot(lock, args.sha, args.ref, args.job_id), sys.stdout, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
