import glob
import os
import re
import numpy as np
import cv2


def _parse_exposure_from_name(fname_no_ext):
    """Parse exposure from common filename conventions.

    Supported examples:
    - "1_4"        -> 1/4 = 0.25
    - "32_1"       -> 32/1 = 32.0
    - "memorial0061_32_1" -> 32/1 = 32.0 (last two numeric parts)
    - "s7708_01.0" -> 1.0   (fallback to last token as float)
    """
    parts = fname_no_ext.split("_")

    # Check for last two parts being integers (handles memorial0061_32_1 format)
    if (
        len(parts) >= 2
        and re.fullmatch(r"\d+", parts[-2])
        and re.fullmatch(r"\d+", parts[-1])
    ):
        num_ = int(parts[-2])
        den_ = int(parts[-1])
        if den_ != 0:
            return num_ / den_

    # Fallback convention: last token is a float exposure time.
    try:
        return float(parts[-1])
    except ValueError:
        return None


def _extract_exposure_from_exif(path):
    """Best-effort EXIF exposure time read.

    Returns a float seconds value, or None if unavailable.
    Uses Pillow only if present; no hard dependency is introduced.
    """
    try:
        from PIL import Image, ExifTags
    except Exception:
        return None

    try:
        img = Image.open(path)
        exif = img.getexif()
        if exif is None:
            return None

        exposure_tag = None
        for tag_id, tag_name in ExifTags.TAGS.items():
            if tag_name == "ExposureTime":
                exposure_tag = tag_id
                break
        if exposure_tag is None:
            return None

        value = exif.get(exposure_tag)
        if value is None:
            return None

        # EXIF exposure is commonly a rational tuple-like value.
        if hasattr(value, "numerator") and hasattr(value, "denominator"):
            if value.denominator == 0:
                return None
            return float(value.numerator) / float(value.denominator)
        if isinstance(value, tuple) and len(value) == 2 and value[1] != 0:
            return float(value[0]) / float(value[1])
        return float(value)
    except Exception:
        return None


def load_images(image_dir, image_ext, root_dir):
    pattern = os.path.join(root_dir, image_dir, image_ext)
    iter_items = sorted(glob.iglob(pattern))

    file_paths = []
    images = []
    exposure_times = []
    for item in iter_items:
        img_bgr = cv2.imread(item)
        if img_bgr is None:
            continue
        file_paths.append(item)
        images.append(cv2.cvtColor(img_bgr, code=cv2.COLOR_BGR2RGB))
        fname = os.path.basename(item)
        exposure_times.append(_parse_exposure_from_name(os.path.splitext(fname)[0]))

    if len(images) == 0:
        raise ValueError(f"No images found for pattern: {pattern}")

    # Try EXIF exposure times for images where filename parsing failed.
    if any(t is None for t in exposure_times):
        exif_times = [_extract_exposure_from_exif(path) for path in file_paths]
        exposure_times = [
            t if t is not None else exif_t
            for t, exif_t in zip(exposure_times, exif_times)
        ]

    # Final fallback: synthetic monotonic exposure times (1..N).
    if any(t is None for t in exposure_times):
        exposure_times = [float(i + 1) for i in range(len(images))]
        print(
            "[load_images] Warning: Could not parse exposure times from filename/EXIF. "
            "Using synthetic exposures 1..N based on sorted filename order."
        )

    images = [
        img
        for _, img in sorted(
            zip(exposure_times, images), key=lambda pair: pair[0], reverse=True
        )
    ]
    B = np.log(sorted(exposure_times, reverse=True))

    return [images, B]
