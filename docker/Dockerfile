FROM brentp/musl-hts-nim:latest AS build


RUN cd / &&    \
    git clone -b master --depth 5 https://github.com/brentp/somalier.git && \
    cd somalier && \
    nimble install -y nimble && \
    /root/.nimble/bin/nimble install -d -y

RUN cd /somalier &&  \
    nim c -d:danger -d:nsb_static -d:release -d:openmp -d:blas=openblas -d:lapack=openblas -o:/usr/bin/somalier src/somalier

FROM alpine:3.21

# for nextflow
RUN apk add --no-cache bash procps

COPY --from=build /usr/bin/somalier /usr/bin/somalier
COPY --from=build /somalier/scripts/ancestry-labels-1kg.tsv /

ENV somalier_ancestry_labels=/ancestry_labels-1kg.tsv


