#!/usr/bin/env python3
import sys

from SigProfilerMatrixGenerator.controllers import cli_controller

if __name__ == "__main__":
    controller = cli_controller.CliController()
    controller.dispatch_install(sys.argv[1:])
