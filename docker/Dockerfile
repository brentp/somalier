FROM brentp/musl-hts-nim:latest AS build


ARG nim_version=2.2.4
ENV PATH=:/root/.nimble/bin:/nim-${nim_version}/bin/:$PATH


RUN cd / && \
    wget -q https://nim-lang.org/download/nim-${nim_version}-linux_x64.tar.xz && \
    tar xf nim-${nim_version}-linux_x64.tar.xz && \
    echo 'PATH=/nim-${nim_version}/bin:$PATH' >> ~/.bashrc && \
    echo 'PATH=/nim-${nim_version}/bin:$PATH' >> ~/.bash_profile && \
    echo 'PATH=/nim-${nim_version}/bin:$PATH' >> /etc/environment  && \
    rm -f nim-${nim_version}-linux_x64.tar.xz


RUN cd / &&    \
    git clone -b master --depth 5 https://github.com/brentp/somalier.git && \
    cd somalier && \
    nimble install -y nimble && \
    nimble install -y arraymancer@0.7.32 && \
    /root/.nimble/bin/nimble install -d -y

# nsb writes to /load assuming it was in docker run -v $(pwd):/load
RUN cd /somalier && mkdir -p /load && \
    nsb -n somalier.nimble -s src/somalier.nim  -- -d:danger -d:nsb_static -d:release -d:blas=openblas -d:lapack=openblas

FROM alpine:3.21

# for nextflow and arraymancer/OpenBLAS runtime loading
RUN apk add --no-cache bash procps openblas

# Copy HTS libraries from build stage
COPY --from=build /usr/local/lib /usr/local/lib
COPY --from=build /usr/local/lib64 /usr/local/lib64
COPY --from=build /usr/lib /usr/lib
COPY --from=build /lib/libcrypto.so.1.1 /lib/libcrypto.so.1.1
COPY --from=build /lib/libssl.so.1.1 /lib/libssl.so.1.1

# Update dynamic linker cache
RUN ldconfig /usr/local/lib /usr/local/lib64 /usr/lib || true
ENV LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:/usr/lib

COPY --from=build /load/somalier /usr/bin/somalier
COPY --from=build /somalier/scripts/ancestry-labels-1kg.tsv /

ENV somalier_ancestry_labels=/ancestry_labels-1kg.tsv
