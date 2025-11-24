"""Replace the links to other crates in documentation strings.

This can be run before building the docs with the `--no-deps` flag to make the links work
in the generated docs
"""

import os
import re


def _insert_doc_links(content):
    content = re.sub(r"(\s)\[ndelement\](\s|\n)", r"\1[ndelement](https://bempp.github.io/ndelement/rust/ndelement/)\1", content)
    return content


def insert_doc_links(folder):
    for file in os.listdir(folder):
        if not file.startswith("."):
            file_path = os.path.join(folder, file)
            if os.path.isdir(file_path):
                insert_doc_links(file_path)
            elif os.path.isfile(file_path) and file.endswith(".rs"):
                with open(file_path) as f:
                    content = f.read()
                with open(file_path, "w") as f:
                    f.write(_insert_doc_links(content))


root_dir = os.path.join(os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    ".."), "..")
insert_doc_links(os.path.join(root_dir, "src"))
insert_doc_links(os.path.join(root_dir, "examples"))
insert_doc_links(os.path.join(root_dir, "tests"))
