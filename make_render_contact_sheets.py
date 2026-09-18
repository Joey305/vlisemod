"""Create temporary contact sheets for visual layout review of rendered pages."""

from pathlib import Path
import sys

from PIL import Image, ImageDraw


def sort_key(path):
    return int(path.stem.split("-")[1])


def main(source_dir, output_prefix):
    source = Path(source_dir)
    pages = sorted(source.glob("page-*.png"), key=sort_key)
    columns, rows = 2, 3
    for sheet_number, start in enumerate(range(0, len(pages), columns * rows), 1):
        chunk = pages[start : start + columns * rows]
        thumbs = []
        for page in chunk:
            image = Image.open(page).convert("RGB")
            image.thumbnail((680, 880))
            thumbs.append((page, image.copy()))
        canvas = Image.new("RGB", (columns * 700, rows * 915), "#d9d9d9")
        draw = ImageDraw.Draw(canvas)
        for index, (page, image) in enumerate(thumbs):
            x = (index % columns) * 700 + 10
            y = (index // columns) * 915 + 25
            canvas.paste(image, (x, y))
            draw.text((x, y - 18), page.stem, fill="black")
        canvas.save(f"{output_prefix}-{sheet_number}.png")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
