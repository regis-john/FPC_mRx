import os
import ast

# Folders you want to scan
FOLDERS = [
    "src/analysis",
    "src/plotting",
    "src/utils",
]

def extract_functions_section(docstring):
    """Return list of (name, description) from the 'Functions:' section."""
    if not docstring:
        return []

    lines = docstring.splitlines()
    functions = []
    in_section = False

    for line in lines:
        stripped = line.strip()

        if stripped.lower().startswith("functions:"):
            in_section = True
            continue

        if in_section:
            if stripped.startswith("-"):
                # Format: "- name: description"
                item = stripped[1:].strip()
                if ":" in item:
                    name, desc = item.split(":", 1)
                    functions.append((name.strip(), desc.strip()))
            else:
                # Stop when section ends
                if stripped:
                    continue
                break

    return functions


def scan_folder(folder):
    """Scan all .py files in a folder and extract function summaries."""
    results = []

    for fname in sorted(os.listdir(folder)):
        if fname.endswith(".py"):
            path = os.path.join(folder, fname)

            with open(path, "r") as f:
                try:
                    module = ast.parse(f.read())
                    docstring = ast.get_docstring(module)
                except SyntaxError:
                    continue

            funcs = extract_functions_section(docstring)
            if funcs:
                results.append((fname, funcs))

    return results


def main():
    for folder in FOLDERS:
        print(f"## {folder}")
        entries = scan_folder(folder)

        if not entries:
            print("(no documented functions)\n")
            continue

        for fname, funcs in entries:
            for name, desc in funcs:
                print(f"- {name}: {desc}")
        print()  # blank line between folders


if __name__ == "__main__":
    main()
