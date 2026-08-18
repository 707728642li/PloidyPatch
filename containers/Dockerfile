# syntax=docker/dockerfile:1.7

ARG BASE_IMAGE=python:3.11.15-slim-bookworm@sha256:d29f48a31a8b408ed19272ca1e7b10ebae13b240a27e862d3d4217c528e2e0c3

FROM ${BASE_IMAGE} AS wheel-builder

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /build
COPY pyproject.toml README.md LICENSE ./
COPY src ./src

# The digest-pinned Python image already contains setuptools 79.0.1. Building
# without isolation and without network access prevents an undeclared build
# dependency from being downloaded into the release candidate.
RUN --network=none python -m pip wheel \
    --no-build-isolation \
    --no-deps \
    --wheel-dir /wheel \
    .

FROM ${BASE_IMAGE} AS runtime

LABEL org.opencontainers.image.title="PloidyPatch core" \
      org.opencontainers.image.description="Review-oriented plant gene-model patch CLI" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.source="https://github.com/707728642li/PloidyPatch"

ENV HOME=/nonexistent \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN groupadd --system --gid 10001 ploidypatch \
    && useradd --system --uid 10001 --gid 10001 \
       --home-dir /nonexistent --shell /usr/sbin/nologin ploidypatch \
    && mkdir -p /work/examples \
    && chown 10001:10001 /work

COPY --from=wheel-builder /wheel/ploidypatch-1.0.0-py3-none-any.whl /tmp/
RUN --network=none python -m pip install \
    --no-deps \
    /tmp/ploidypatch-1.0.0-py3-none-any.whl \
    && rm /tmp/ploidypatch-1.0.0-py3-none-any.whl

USER 10001:10001
WORKDIR /work

ENTRYPOINT ["ploidypatch"]
CMD ["--help"]
