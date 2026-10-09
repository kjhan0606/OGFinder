#!/bin/bash
# Stage a system OpenSSL 3 for the ds9 link.
# The vendored openssl/ tree is 1.0.2u and is not configured or compiled.
# Destination layout matches tls 1.6 configure: $dest/include/openssl and $dest/lib.
set -euo pipefail

dest="${1:?usage: stage_openssl3.sh DEST}"

pick_lib() {
    local root="$1"
    local d
    for d in "$root/lib64" "$root/lib" "$root"; do
        if [ -e "$d/libssl.so" ] || [ -e "$d/libssl.so.3" ] \
            || [ -e "$d/libssl.dylib" ] || [ -e "$d/libssl.3.dylib" ] \
            || [ -e "$d/libssl-3.dll" ] || [ -e "$d/libssl.dll" ]; then
            printf '%s\n' "$d"
            return 0
        fi
    done
    return 1
}

refuse_old() {
    local root="$1"
    if [ -e "$root/lib/ssleay32.dll" ] || [ -e "$root/lib/libeay32.dll" ] \
        || [ -e "$root/ssleay32.dll" ] || [ -e "$root/libeay32.dll" ]; then
        if ! pick_lib "$root" >/dev/null 2>&1; then
            echo "OpenSSL 1.0 (ssleay32/libeay32) is not used. OpenSSL 3 is required." >&2
            exit 1
        fi
    fi
}

inc=""
lib=""
if [ -n "${OPENSSL_DIR:-}" ]; then
    refuse_old "$OPENSSL_DIR"
    if [ ! -f "$OPENSSL_DIR/include/openssl/opensslv.h" ]; then
        echo "OPENSSL_DIR has no include/openssl/opensslv.h: $OPENSSL_DIR" >&2
        exit 1
    fi
    inc="$OPENSSL_DIR/include"
    lib="$(pick_lib "$OPENSSL_DIR")" || {
        echo "OPENSSL_DIR has no libssl: $OPENSSL_DIR" >&2
        exit 1
    }
elif [ -f /usr/include/openssl/opensslv.h ]; then
    inc=/usr/include
    lib="$(pick_lib /usr)" || {
        echo "Found /usr/include/openssl but no libssl under /usr/lib64 or /usr/lib." >&2
        exit 1
    }
else
    echo "No system OpenSSL headers. Set OPENSSL_DIR to an OpenSSL 3 prefix." >&2
    echo "The vendored openssl/ 1.0.2u tree is not built." >&2
    exit 1
fi

major="$(sed -n 's/^#[[:space:]]*define[[:space:]]\{1,\}OPENSSL_VERSION_MAJOR[[:space:]]\{1,\}\([0-9][0-9]*\).*/\1/p' \
    "$inc/openssl/opensslv.h" | head -n 1)"
case "$major" in
    ''|*[!0-9]*)
        echo "OpenSSL 3 is required. $inc/openssl/opensslv.h has no OPENSSL_VERSION_MAJOR." >&2
        echo "The vendored openssl/ 1.0.2u tree is not built." >&2
        exit 1
        ;;
esac
if [ "$major" -lt 3 ]; then
    echo "OpenSSL $major is not used. OpenSSL 3 is required." >&2
    echo "The vendored openssl/ 1.0.2u tree is not built." >&2
    exit 1
fi

inc="$(readlink -f "$inc")"
lib="$(readlink -f "$lib")"
rm -rf "$dest"
mkdir -p "$dest"
ln -s "$inc" "$dest/include"
ln -s "$lib" "$dest/lib"
echo "staged OpenSSL $major include=$inc lib=$lib -> $dest"
