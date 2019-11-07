FROM brentp/musl-hts-nim:latest

RUN cd / &&   \
    git clone -b master --depth 1 git://github.com/brentp/somalier.git && \
    cd somalier && \
    nimble install -d -y && \
    nim c -d:danger -o:/usr/bin/somalier -d:openmp -d:openblas --passC:-flto src/somalier && \
    rm -rf /somalier && somalier --help

