#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 IMAGE EXAMPLE_DIR OUTPUT_DIR" >&2
  exit 2
fi

image=$1
example_dir=$2
output_dir=$3
working_dir="${output_dir}.working"

for executable in docker python3 sha256sum; do
  command -v "$executable" >/dev/null 2>&1 || {
    echo "missing executable: $executable" >&2
    exit 1
  }
done

[[ -d "$example_dir" && ! -L "$example_dir" ]] || {
  echo "example directory must be a regular directory: $example_dir" >&2
  exit 1
}
example_dir=$(cd "$example_dir" && pwd -P)

if [[ -e "$output_dir" || -e "$working_dir" ]]; then
  echo "refusing to overwrite container smoke output: $output_dir" >&2
  exit 1
fi
mkdir -p "$(dirname "$output_dir")"
mkdir "$working_dir"

docker image inspect "$image" > "$working_dir/image_inspect.json"
docker run --rm --read-only --network none \
  "$image" --version > "$working_dir/version.txt"
docker run --rm --read-only --network none \
  "$image" --help > "$working_dir/cli_help.txt"

for family in audit baseline benchmark evidence graph normalize patch; do
  grep -qw "$family" "$working_dir/cli_help.txt" || {
    echo "container help is missing command family: $family" >&2
    exit 1
  }
done

docker run --rm \
  --read-only \
  --network none \
  --tmpfs /tmp:rw,noexec,nosuid,size=16m \
  --mount "type=bind,src=${example_dir},dst=/work/examples/minimal_reviewed_patch,readonly" \
  --entrypoint python \
  "$image" \
  /work/examples/minimal_reviewed_patch/run_example.py \
  --input-dir /work/examples/minimal_reviewed_patch \
  --output-dir /tmp/reviewed-patch \
  > "$working_dir/example_stdout.json"

python3 - \
  "$image" \
  "$working_dir/image_inspect.json" \
  "$working_dir/version.txt" \
  "$working_dir/example_stdout.json" \
  "$working_dir/container_smoke.json" <<'PY'
import json
import sys
from pathlib import Path

image_reference, inspect_path, version_path, example_path, output_path = sys.argv[1:]
inspection = json.loads(Path(inspect_path).read_text(encoding="utf-8"))
if not isinstance(inspection, list) or len(inspection) != 1:
    raise ValueError("docker inspect must return exactly one image")
image = inspection[0]
config = image.get("Config") or {}
user = config.get("User")
entrypoint = config.get("Entrypoint")
version = Path(version_path).read_text(encoding="utf-8").strip()
example = json.loads(Path(example_path).read_text(encoding="utf-8"))

checks = {
    "runtime_user_is_10001": user == "10001:10001",
    "entrypoint_is_ploidypatch": entrypoint == ["ploidypatch"],
    "version_is_1_1_0": version == "1.0.1",
    "accepted_additions_are_two": example.get("accepted_additions") == 2,
    "automatic_approval_is_false": example.get("automatic_approval") is False,
    "byte_identical_reversion": example.get("byte_identical_reversion") is True,
    "source_reverted_sha_match": (
        example.get("source_sha256") == example.get("reverted_sha256")
    ),
}
failed = sorted(name for name, passed in checks.items() if not passed)
if failed:
    raise ValueError("container smoke checks failed: " + ", ".join(failed))

report = {
    "schema_version": "ploidypatch.container_smoke.v1",
    "image_reference": image_reference,
    "image_id": image.get("Id"),
    "image_size_bytes": image.get("Size"),
    "architecture": image.get("Architecture"),
    "os": image.get("Os"),
    "runtime_user": user,
    "entrypoint": entrypoint,
    "version": version,
    "runtime_read_only": True,
    "runtime_network": "none",
    "checks": checks,
    "example": example,
}
Path(output_path).write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
    newline="",
)
PY

(
  cd "$working_dir"
  find . -type f ! -name SHA256SUMS -printf '%P\0' \
    | sort -z \
    | xargs -0 sha256sum
) > "$working_dir/SHA256SUMS"

mv "$working_dir" "$output_dir"
cat "$output_dir/container_smoke.json"
