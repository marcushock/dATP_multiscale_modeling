FROM nvidia/cuda:12.2.0-devel-ubuntu22.04
RUN apt update && apt install -y libboost-system-dev libboost-thread-dev
