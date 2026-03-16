"""Download and cache pre-built translation bundles from GitHub Releases."""

from __future__ import annotations

import json
import logging
import urllib.request
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)

GITHUB_REPO = "maltesie/gtdb-translate"
GITHUB_API_RELEASES = (
    f"https://api.github.com/repos/{GITHUB_REPO}/releases"
)

DEFAULT_CACHE_DIR = Path.home() / ".cache" / "gtdb_translate"
BUNDLE_FILENAME_TEMPLATE = "gtdb_translate_{version}.msgpack.zst"


def _resolve_latest_version() -> str:
    """Query the GitHub API for the latest release that has a bundle asset.

    Iterates through releases (newest first) and returns the tag of the
    first one whose assets include a file matching
    ``gtdb_translate_<tag>.msgpack.zst``.

    Returns
    -------
    str
        The tag name of the latest bundle release (e.g. ``"r226"``).

    Raises
    ------
    RuntimeError
        If the API request fails or no release with a bundle is found.
    """
    req = urllib.request.Request(
        GITHUB_API_RELEASES,
        headers={"Accept": "application/vnd.github.v3+json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            releases = json.loads(resp.read().decode())
    except Exception as exc:
        raise RuntimeError(
            f"Failed to query releases from {GITHUB_API_RELEASES}. "
            f"You can specify a version explicitly instead.\n\n"
            f"Original error: {exc}"
        ) from exc

    for release in releases:
        tag = release.get("tag_name", "")
        expected_asset = BUNDLE_FILENAME_TEMPLATE.format(version=tag)
        asset_names = [a.get("name", "") for a in release.get("assets", [])]
        if expected_asset in asset_names:
            logger.info("Latest bundle release: %s", tag)
            return tag

    raise RuntimeError(
        f"No release with a translation bundle found in "
        f"https://github.com/{GITHUB_REPO}/releases. "
        f"You can specify a version explicitly instead."
    )


def _asset_download_url(version: str) -> str:
    """Construct the direct download URL for a release asset."""
    filename = BUNDLE_FILENAME_TEMPLATE.format(version=version)
    return (
        f"https://github.com/{GITHUB_REPO}/releases/download/"
        f"{version}/{filename}"
    )


def ensure_bundle(
    version: Optional[str] = None,
    cache_dir: Optional[Union[str, Path]] = None,
    force: bool = False,
) -> Path:
    """Return the local path to a translation bundle, downloading if needed.

    Parameters
    ----------
    version : str, optional
        GTDB release version (e.g. ``"r226"``).  If ``None``, the latest
        release is resolved automatically via the GitHub API.
    cache_dir : str or Path, optional
        Override the default cache directory (``~/.cache/gtdb_translate/``).
    force : bool
        If ``True``, re-download even if the file already exists locally.

    Returns
    -------
    Path
        Path to the cached ``.msgpack.zst`` file.
    """
    if version is None:
        version = _resolve_latest_version()

    cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)
    filename = BUNDLE_FILENAME_TEMPLATE.format(version=version)
    local_path = cache_dir / filename

    if local_path.exists() and not force:
        logger.debug("Using cached bundle: %s", local_path)
        return local_path

    url = _asset_download_url(version)
    logger.info("Downloading bundle from %s ...", url)

    try:
        urllib.request.urlretrieve(url, local_path)
    except Exception as exc:
        if local_path.exists():
            local_path.unlink()
        raise RuntimeError(
            f"Failed to download translation bundle for {version} from "
            f"{url}. Make sure the release exists and you have internet "
            f"access.\n\nOriginal error: {exc}"
        ) from exc

    logger.info("Saved to %s (%s)", local_path, version)
    return local_path
