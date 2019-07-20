FROM brentp/musl-hts-nim:latest

RUN cd / &&      \
    git clone -b master --depth 1 git://github.com/brentp/somalier.git && \
    cd somalier && \
    nimble install -y && \
    nim c -d:danger -o:/usr/bin/somalier --passC:-flto src/somalier && \
    rm -rf /somalier && somalier --help

