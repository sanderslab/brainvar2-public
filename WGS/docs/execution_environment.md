# Execution environment

This page separates NVIDIA's version-matched minimum requirements from the
resource allocation used successfully for BrainVar2 WGS. The registry and
local execution checks were repeated on 2026-07-15.

## Exact Parabricks image

The exact production image remains downloadable from NVIDIA NGC:

```text
nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1
```

The workflow pins its immutable manifest digest:

```text
nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1@sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e
```

An anonymous remote manifest inspection resolved the tag to that digest. HTTP
HEAD checks reached the manifest, image configuration, and every image layer.
The image is a Docker schema-v2 `linux/amd64` image, and the
[NVIDIA NGC catalog](https://catalog.ngc.nvidia.com/orgs/nvidia/clara/containers/clara-parabricks/4.0.1-1/tags)
lists the `4.0.1-1` release.

NVIDIA documents that guest users can download public NGC catalog entities
without signing in. If a site requires authenticated NGC access, use an NGC
API key as documented in the
[NGC Catalog access guide](https://docs.nvidia.com/ngc/latest/ngc-catalog-user-guide.html#accessing-ngc-software).
NVIDIA's
[4.0 licensing note](https://docs.nvidia.com/clara/parabricks/about-parabricks/release-notes/4-0-0-1-release-notes)
also states that a separate Parabricks license is no longer required starting
with version 4.0.0-1; enterprise support is optional.

The archived 4.0.1 guide shows a `v4.0.1-1` tag. That prefixed tag no longer
resolves, while the unprefixed `4.0.1-1` tag above does. Use
the workflow's pinned unprefixed reference:

```bash
IMAGE='nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1@sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e'
docker manifest inspect "$IMAGE"
docker pull "$IMAGE"
docker run --rm --gpus all "$IMAGE" pbrun version
```

The final command completed locally and reported `pbrun: 4.0.1-1`. An
unofficial `jlalli/parabricks` duplicate is not needed; the exact production
image is available directly from NVIDIA.

## Version-matched NVIDIA requirements

The archived
[Parabricks 4.0.1 installation requirements](https://docs.nvidia.com/clara/parabricks/4.0.1/gettingstarted.html)
specify:

| Component | NVIDIA Parabricks 4.0.1 requirement |
| --- | --- |
| GPU | At least 16 GB GPU RAM; the guide lists CUDA architectures 60, 70, 75, and 80 and reports testing on V100, T4, A10, A30, A40, A100, and A6000 GPUs |
| Two-GPU host | At least 100 GB CPU RAM and 24 CPU threads |
| Four-GPU host | At least 196 GB CPU RAM and 32 CPU threads |
| Eight-GPU host | At least 392 GB CPU RAM and 48 CPU threads |
| Driver | NVIDIA driver later than 465.32.x |
| Runtime | Linux supporting `nvidia-docker2`, with Docker 20.10 or later |
| Virtualization | vGPU and MIG are not supported by this release |

The same guide recommends fast local SSD storage for inputs, outputs, and
temporary files. These are release-specific requirements; requirements in the
current Parabricks documentation apply to newer images and are not substituted
for the archived 4.0.1 contract.

## BrainVar2-tested allocation

The retained production run record says each per-library job requested two
GPUs, 24 CPU cores, 200 GB host RAM, and 400 GB disk, and notes that some
libraries needed additional RAM or disk. A later complete three-FASTQ-chunk
BrainVar2 sample reproduced the retained gVCF and index byte-for-byte using:

- two NVIDIA L4 GPUs with 23,034 MiB each;
- NVIDIA driver 580.105.08;
- Docker 29.1.3 with NVIDIA Container Toolkit 1.18.1;
- 24 CPU threads and a 200 GB Nextflow memory allocation.

In that run, `fq2bam` completed in 1h02m05s with 121.9 GB peak host RSS, and
HaplotypeCaller completed in 40m09s with 109 GB peak host RSS. The L4 is newer
than the GPU list in NVIDIA's archived 4.0.1 guide, so this is empirical
compatibility evidence rather than a claim of version-matched vendor support.

The workflow therefore requests a compatible NVIDIA GPU resource without
requiring the two-GPU allocation used by the equivalence run. Sites may assign
additional GPUs through their local Nextflow configuration. The retained run
used 24 CPU threads, 200 GB host RAM, and 400 GB disk; increase RAM and scratch
for unusually large or multi-chunk libraries. The joint-genotyping/VQSR stage
is CPU-only and uses the separate GATK container; it does not require a GPU.

## Preflight checks

Run these before protected data are staged:

```bash
nvidia-smi
docker version
docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi

IMAGE='nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1@sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e'
docker pull "$IMAGE"
docker run --rm --gpus all "$IMAGE" pbrun version
```

Run the preflight and synthetic stub workflow once in each new execution
environment before staging controlled BrainVar2 FASTQs.
