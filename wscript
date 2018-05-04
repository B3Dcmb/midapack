#! /usr/bin/env python
# encoding: utf-8


def build(ctx):

    lstlib = ['M', 'FFTW3']
    if ctx.env.FFTWMT:
        lstlib += ['FFTW3MULTITHREAD']
    if ctx.env.OMP:
        lstlib += ['OMP', 'PTHREAD']

    ## Toeplitz core library
    ctx.new_task_gen(features = 'c cstlib',
                     source = ['./toeplitz.c','./toeplitz_seq.c','./toeplitz_nofft.c','./toeplitz_gappy.c','./toeplitz_block.c', 'toeplitz_devtools.c', 'toeplitz_rshp.c', 'toeplitz_params.c', 'toeplitz_utils.c', 'toeplitz_wizard.c'],
                     target='toeplitz', 
                     use = lstlib,
		     cflags = ctx.env['CFLAGS'])

