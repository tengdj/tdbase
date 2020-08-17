/* Based on the following file, stripped down and changed names by MS
   Original header follows. Contrary to below please send bug reports
   to michael@compressconsult.com
   The bytecount field is misused to store the number of free bits.
   To activate this code place the file in the same directory as the
   rangecod.c and #define RENORM95 in rangecod.c
   According to the notice below this file is provided for research
   purposes on http://www.compressconsult.com/rangecoder/
*/

/******************************************************************************
File:           arith.c

Authors:        John Carpinelli   (johnfc@ecr.mu.oz.au)
                Wayne Salamonsen  (wbs@mundil.cs.mu.oz.au)
                Lang Stuiver      (langs@cs.mu.oz.au)
                Radford Neal      (radford@ai.toronto.edu)

Purpose:        Data compression using a revised arithmetic coding method.

Based on:       A. Moffat, R. Neal, I.H. Witten, "Arithmetic Coding Revisted",
                Proc. IEEE Data Compression Conference, Snowbird, Utah,
                March 1995.

                Low-Precision Arithmetic Coding Implementation by
                Radford M. Neal



Copyright 1995 John Carpinelli and Wayne Salamonsen, All Rights Reserved.
Copyright 1996 Lang Stuiver, All Rights Reserved.

These programs are supplied free of charge for research purposes only,
and may not sold or incorporated into any commercial product.  There is
ABSOLUTELY NO WARRANTY of any sort, nor any undertaking that they are
fit for ANY PURPOSE WHATSOEVER.  Use them at your own risk.  If you do
happen to find a bug, or have modifications to suggest, please report
the same to Alistair Moffat, alistair@cs.mu.oz.au.  The copyright
notice above and this statement of conditions must remain an integral
part of each and every copy made of these files.

******************************************************************************/

char coderversion[] = "Range Coder 1.1 with \"DCC95 renormalisation, Moffat et al.\". Modified by MS";

#define         Half            (Top_value>>1)
#define         Quarter         (Half>>1)

#define OUTPUT_BIT(rco,b)                \
do {if (!(rco)->bytecount--)             \
    {   outbyte((rco),(rco)->buffer);    \
        (rco)->bytecount = 7;            \
    }                                    \
    (rco)->buffer = (rco)->buffer<<1 | b;\
} while (0)


#define INPUT_BIT(rco,target)            \
do {if (!(rco)->bytecount)               \
    {   (rco)->buffer = inbyte(rco);     \
        (rco)->bytecount = 8;            \
    }                                    \
    (rco)->bytecount--;                  \
    (target) = (target)<<1 | (((rco)->buffer>>(rco)->bytecount)&1);\
} while (0)


/*
 * BIT_PLUS_FOLLOW(b)
 * responsible for outputting the bit passed to it and an opposite number of
 * bit equal to the value stored in help (bits_outstanding)
 */
#define BIT_PLUS_FOLLOW(rco,b)           \
do {OUTPUT_BIT((rco),(b));               \
    while ((rco)->help > 0)              \
    {                                    \
        OUTPUT_BIT((rco),!(b));          \
        (rco)->help--;                   \
    }                                    \
} while (0)


/*
 * ENCODE_RENORMALISE
 * output code bits until the range has been expanded
 * to above QUARTER
 */
static Inline void enc_normalize( rangecoder *rc )
{   while (rc->range <= Quarter)
    {
        if (rc->low >= Half)
        {
            BIT_PLUS_FOLLOW(rc,1);
            rc->low -= Half;
        }
        else if (rc->low+rc->range <= Half)
        {
            BIT_PLUS_FOLLOW(rc,0);
        }
        else
        {
            rc->help++;
            rc->low -= Quarter;
        }
        rc->low <<= 1;
        rc->range <<= 1;
    }
}


uint4 done_encoding( rangecoder *rc )
{   int i,b;
    rc->low += rc->range>>1;
    enc_normalize( rc );
    for (i = CODE_BITS-1; i>0; i--)
    {   b = (rc->low>>i)&1;
        BIT_PLUS_FOLLOW(rc,b);
    }
    b = rc->low & 1;
    BIT_PLUS_FOLLOW(rc,b);
    outbyte(rc,rc->buffer<<rc->bytecount);
    return rc->bytecount;
}


/* Start the decoder                                         */
/* rc is the range coder to be used                          */
/* returns the char from start_encoding or EOF               */
int start_decoding( rangecoder *rc )
{   int c = inbyte(rc);
    if (c==EOF)
        return EOF;
    rc->bytecount = 0;
    rc->low = 0;
    rc->range = 1;
    return c;
}

/*
 * DECODE_RENORMALISE
 * input code bits until range has been expanded to
 * more than QUARTER. Mimics encoder.
 */
static Inline void dec_normalize( rangecoder *rc )
{   while (rc->range <= Quarter)
    {
        rc->range <<= 1;
        INPUT_BIT(rc,rc->low);
    }
}
