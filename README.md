# `gs.c` = A single C file Gaussian Splatting renderer

This project is an accurate(to the best of my knowledge) implementation of the [Gaussian Splatting](https://repo-sam.inria.fr/fungraph/3d-gaussian-splatting/) project's rendering pipeline. I wrote this program to learn both `C` and `GS` from scratch. Inspired by [Karpathy's llm.c](https://github.com/karpathy/llm.c), I wrote it in a single C file with no external dependency. Just use the `-lm` math flag while compiling. 


## Requirements
- C compiler (`gcc`)
- ImageMagick

## How to run
- You can use the `run.sh` script. You can also manually compile and run it. 
```bash
chmod +x run.sh # might be required
./run.sh
```
- The final render is the file `render_gs_c.png`
- You can edit the following line in `gs.c` to render different views in the dataset
```c
// change 10 to different view. This num is the image ID in the cameras.json or the dataset
#define TEST_IMG_ID 10 
```
## References
- [Official Implementation by the GS authors](https://github.com/graphdeco-inria/gaussian-splatting/tree/main)
- [hbb1/torch-splatting](https://github.com/hbb1/torch-splatting)
- [MrNeRF's C++/CUDA implementation](https://github.com/MrNeRF/gaussian-splatting-cuda)
- [NeRF studio's gsplat](https://github.com/nerfstudio-project/gsplat/tree/bbc9e98fa9add8b355b91b69096ec3b1271a9f46)
- [you will never ask about pointers again after watching this video](https://www.youtube.com/watch?v=2ybLD6_2gKM)
- [OpenSplat](https://github.com/pierotofy/OpenSplat/blob/main/rasterizer/gsplat-cpu/gsplat_cpu.cpp)

## TODO
- [ ] Port to pthreads or openMP for multithreading
- [ ] Remove debugging code

## DISCLAIMER/WARRANTY
I give no warranty for this repo. Be precausious and use it. It might break your computer and Im not responsible for it. This repo is meant for research purposes only. Im a beginner C programmer, so this repo might not be the best memory-safe efficient code. Im still learning. I humbly accept any mistakes and stupidity in the code. Im open to pull requests and suggestions. I will take time to implement/merge as Im also a full time student and I need to understand your suggestions completely before adding them.