/*
 * Copyright (C) 1997-2000 Matt Newman <matt@novadigm.com>
 *
 * $Header: /cvsroot/tls/tls/tlsBIO.c,v 1.8 2004/03/24 05:22:53 razzell Exp $
 *
 * Provides BIO layer to interface openssl to Tcl.
 *
 * OpenSSL 1.1 and 3 keep BIO opaque. The product links OpenSSL 3, so this
 * file uses BIO_meth_* and BIO_get_data instead of the OpenSSL 1.0 fields.
 */

#include "tlsInt.h"

/*
 * Forward declarations
 */

static int BioWrite(BIO *h, const char *buf, int num);
static int BioRead(BIO *h, char *buf, int num);
static int BioPuts(BIO *h, const char *str);
static long BioCtrl(BIO *h, int cmd, long arg1, void *ptr);
static int BioNew(BIO *h);
static int BioFree(BIO *h);

static BIO_METHOD *BioMethods = NULL;

static BIO_METHOD *BioMethod(void)
{
    if (BioMethods == NULL) {
	BioMethods = BIO_meth_new(BIO_TYPE_TCL, "tcl");
	if (BioMethods == NULL) {
	    return NULL;
	}
	BIO_meth_set_write(BioMethods, BioWrite);
	BIO_meth_set_read(BioMethods, BioRead);
	BIO_meth_set_puts(BioMethods, BioPuts);
	BIO_meth_set_ctrl(BioMethods, BioCtrl);
	BIO_meth_set_create(BioMethods, BioNew);
	BIO_meth_set_destroy(BioMethods, BioFree);
    }
    return BioMethods;
}

BIO *BIO_new_tcl(State *statePtr, int flags)
{
    BIO *bio;
    BIO_METHOD *method = BioMethod();

    if (method == NULL) {
	return NULL;
    }
    bio = BIO_new(method);
    if (bio == NULL) {
	return NULL;
    }
    BIO_set_data(bio, (char *)statePtr);
    BIO_set_init(bio, 1);
    BIO_set_shutdown(bio, flags);

    return bio;
}

BIO_METHOD *BIO_s_tcl()
{
    return BioMethod();
}

static int BioWrite (BIO *bio, const char *buf, int bufLen)
{
    Tcl_Channel chan = Tls_GetParent((State *)BIO_get_data(bio));
    int ret;

    dprintf(stderr,"\nBioWrite(0x%p, <buf>, %d) [0x%p]",
	    (void*)bio, bufLen, (void*) chan);

    if (channelTypeVersion == TLS_CHANNEL_VERSION_2) {
	ret = Tcl_WriteRaw(chan, buf, bufLen);
    } else {
	ret = Tcl_Write(chan, buf, bufLen);
    }

    dprintf(stderr,"\n[0x%p] BioWrite(%d) -> %d [%d.%d]",
	    (void*) chan, bufLen, ret, Tcl_Eof(chan), Tcl_GetErrno());

    BIO_clear_flags(bio, BIO_FLAGS_WRITE|BIO_FLAGS_SHOULD_RETRY);

    if (ret == 0) {
	if (!Tcl_Eof(chan)) {
	    BIO_set_retry_write(bio);
	    ret = -1;
	}
    }
    if (BIO_should_read(bio)) {
	BIO_set_retry_read(bio);
    }
    return ret;
}

static int BioRead (BIO *bio, char *buf, int bufLen)
{
    Tcl_Channel chan = Tls_GetParent((State *)BIO_get_data(bio));
    int ret = 0;

    dprintf(stderr,"\nBioRead(0x%p, <buf>, %d) [0x%p]",
	    (void*) bio, bufLen, (void*) chan);

    if (buf == NULL) return 0;

    if (channelTypeVersion == TLS_CHANNEL_VERSION_2) {
	ret = Tcl_ReadRaw(chan, buf, bufLen);
    } else {
	ret = Tcl_Read(chan, buf, bufLen);
    }

    dprintf(stderr,"\n[0x%p] BioRead(%d) -> %d [%d.%d]",
	    (void*) chan, bufLen, ret, Tcl_Eof(chan), Tcl_GetErrno());

    BIO_clear_flags(bio, BIO_FLAGS_READ|BIO_FLAGS_SHOULD_RETRY);

    if (ret == 0) {
	if (!Tcl_Eof(chan)) {
	    BIO_set_retry_read(bio);
	    ret = -1;
	}
    }
    if (BIO_should_write(bio)) {
	BIO_set_retry_write(bio);
    }
    return ret;
}

static int BioPuts (BIO *bio, const char *str)
{
    return BioWrite(bio, str, (int) strlen(str));
}

static long BioCtrl(BIO *bio, int cmd, long num, void *ptr)
{
    Tcl_Channel chan = Tls_GetParent((State *)BIO_get_data(bio));
    long ret = 1;

    dprintf(stderr,"\nBioCtrl(0x%p, 0x%x, 0x%p, 0x%p)",
	    (void*)bio, (unsigned int)cmd, (void*)num, ptr);

    switch (cmd) {
    case BIO_CTRL_RESET:
	num = 0;
    case BIO_C_FILE_SEEK:
    case BIO_C_FILE_TELL:
	ret = 0;
	break;
    case BIO_CTRL_INFO:
	ret = 1;
	break;
    case BIO_C_SET_FD:
	BioFree(bio);
	BIO_set_data(bio, *((char **)ptr));
	BIO_set_shutdown(bio, (int)num);
	BIO_set_init(bio, 1);
	break;
    case BIO_C_GET_FD:
	if (BIO_get_init(bio)) {
	    if (ptr != NULL) {
		*(int *)ptr = 0;
	    }
	    ret = 0;
	} else {
	    ret = -1;
	}
	break;
    case BIO_CTRL_GET_CLOSE:
	ret = BIO_get_shutdown(bio);
	break;
    case BIO_CTRL_SET_CLOSE:
	BIO_set_shutdown(bio, (int)num);
	break;
    case BIO_CTRL_EOF:
	dprintf(stderr, "BIO_CTRL_EOF\n");
	ret = Tcl_Eof(chan);
	break;
    case BIO_CTRL_PENDING:
	ret = (Tcl_InputBuffered(chan) ? 1 : 0);
	dprintf(stderr, "BIO_CTRL_PENDING(%d)\n", (int) ret);
	break;
    case BIO_CTRL_WPENDING:
	ret = 0;
	break;
    case BIO_CTRL_DUP:
	break;
    case BIO_CTRL_FLUSH:
	dprintf(stderr, "BIO_CTRL_FLUSH\n");
	if (channelTypeVersion == TLS_CHANNEL_VERSION_2) {
	    ret = ((Tcl_WriteRaw(chan, "", 0) >= 0) ? 1 : -1);
	} else {
	    ret = ((Tcl_Flush(chan) == TCL_OK) ? 1 : -1);
	}
	break;
    default:
	ret = 0;
	break;
    }
    return(ret);
}

static int BioNew(BIO *bio)
{
    BIO_set_init(bio, 0);
    BIO_set_data(bio, NULL);
    return 1;
}

static int BioFree(BIO *bio)
{
    if (bio == NULL) {
	return 0;
    }

    if (BIO_get_shutdown(bio)) {
	if (BIO_get_init(bio)) {
	    BIO_set_init(bio, 0);
	}
    }
    return 1;
}
