import tomllib
with open("Manifest.toml", "rb") as f:
    print(tomllib.load(f)["julia_version"])
