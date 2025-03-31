## Deconvolution of Slide-seq Data

Deconvolution Methods:
- RCTD
- STdeconvolve

#### 1. Environment Setup
Run these commands locally then transfer to hpc.

```
docker buildx build -t st_container --platform linux/amr64 -o path/to/dockerfile .
```

```
docker save st_container -o st_container.tar
```

```
scp path/to/st_container.tar username@remote_host:///path/to/destination/
```
