# -*- coding: utf-8 -*-
"""
Created on 2026-08-28T16:23:41-04:00

@author: nate
"""
import argh
from loguru import logger

import lex-earliest-seqs


def main():
    logger.info(__name__)

def cli():
    parser = argh.ArghParser()
    parser.add_commands([
            main
    ])
    parser.dispatch()

    # Only one entrypoint
    #argh.dispatch_command(main)