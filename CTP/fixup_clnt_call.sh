#!/bin/sh

TMPF=tmp_clnt.c
if [ -r "$1" ]; then
    if ! grep -q clnt_call.*xdrproc_t "$1"; then
        sed -E 's|clnt_call\(([^,]+, )([^,]+, )([^,]+, )([^,]+, )([^,]+, )(.*)(\) .*)$|clnt_call(\1\2(xdrproc_t)\3\4(xdrproc_t)\5\6\7|' "$1" > $TMPF
        mv -f $TMPF "$1"
    fi
fi
