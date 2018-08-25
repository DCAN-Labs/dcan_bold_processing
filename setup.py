import argparse


def _cli():
    parser = generate_parser()
    args = parser.parse_args()

    return interface(args)


def generate_parser(parser=None):
    parser = argparse.ArgumentParser(
        prog='setup.py',
        description="""
        compiles matlab code for use."""
    )
    parser.add_argument('matlab_compiler', nargs=1, type=str,
                        help='path to the matlab compiler binary.')

    return parser


def interface():
    pass


if __name__ == '__main__':
    _cli()
