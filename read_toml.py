import tomllib

with open("Manifest.toml", "rb") as f:
    data = tomllib.load(f)
    print(data["julia_version"])
