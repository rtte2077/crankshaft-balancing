import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
READMES = ["README.md", "README.en.md"]
LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)\)")


@pytest.mark.parametrize("name", READMES)
def test_relative_links_and_images_exist(name):
    text = (ROOT / name).read_text(encoding="utf-8")
    missing = []
    for target in LINK.findall(text):
        if target.startswith(("http://", "https://", "#", "mailto:")):
            continue
        path = target.split("#")[0]
        if path and not (ROOT / path).exists():
            missing.append(target)
    assert not missing, f"{name}: нет файлов {missing}"


def test_readmes_link_to_each_other():
    assert "README.en.md" in (ROOT / "README.md").read_text(encoding="utf-8")
    assert "README.md" in (ROOT / "README.en.md").read_text(encoding="utf-8")


@pytest.mark.parametrize("name", READMES)
def test_no_latex_macros_rejected_by_github(name):
    # GitHub (MathJax) не принимает \operatorname — формула превращается в сообщение об ошибке
    text = (ROOT / name).read_text(encoding="utf-8")
    assert "\\operatorname" not in text
