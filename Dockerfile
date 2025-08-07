# Stage 1: Build stage
FROM python:3.13-alpine AS build

WORKDIR /app

# Install build dependencies (including build-base, python3-dev)
RUN apk add --no-cache --virtual .build-deps \
    build-base \
    python3-dev \
    && pip install --upgrade pip setuptools wheel

# Copy only necessary files
COPY pyproject.toml CHANGELOG.md README.md MANIFEST.in LICENSE /app/
COPY enzymm /app/enzymm

# build wheels into /wheels
RUN pip wheel . -w /wheels

# Stage 2: Runtime stage
FROM python:3.13-alpine AS run

ARG VERSION=latest
LABEL org.opencontainers.image.title="enzymm"
LABEL org.opencontainers.image.description="Enzyme Motif Miner"
LABEL org.opencontainers.image.version=$VERSION
LABEL org.opencontainers.image.authors="Raymund Hackett <r.e.hackett@lumc.nl>"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.url="https://github.com/rayhackett/enzymm"
LABEL org.opencontainers.image.source="https://github.com/rayhackett/enzymm"

WORKDIR /app

# Copy built wheels from build stage
COPY --from=build /wheels /wheels

# Install wheels without build tools
RUN pip install --no-cache-dir /wheels/*

# Entry point in pyproject.toml set with [project.scripts]\n enzymm = "enzymm._cli:main"
CMD ["enzymm"]

