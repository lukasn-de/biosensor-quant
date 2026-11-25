"""Command line: `biosensor-quant run`, `biosensor-quant entropy`, `biosensor-quant demo`.

Every setting in config.PARAMS is available as an option (underscores become dashes).
"""

import argparse

from .config import DEFAULTS, PARAMS, PRESETS, build_settings
from .pipeline import entropy_table, run_acquisition


def _bool(text):
    if text.lower() in ("true", "yes", "1", "on"):
        return True
    if text.lower() in ("false", "no", "0", "off"):
        return False
    raise argparse.ArgumentTypeError(f"expected true/false, got {text!r}")


def _size(text):
    h, w = (int(x) for x in text.replace("x", ",").split(","))
    return [h, w]


def _add_setting_options(parser):
    group = parser.add_argument_group("settings (defaults follow the thesis method)")
    for name, (default, kind, help_text) in PARAMS.items():
        kwargs = {"default": None, "help": f"{help_text} [default: {default}]"}
        if isinstance(kind, list):
            kwargs["choices"] = kind
        elif kind is bool:
            kwargs["type"], kwargs["metavar"] = _bool, "true|false"
        elif kind == "size":
            kwargs["type"], kwargs["metavar"] = _size, "H,W"
        else:
            kwargs["type"] = kind
        group.add_argument("--" + name.replace("_", "-"), dest=name, **kwargs)


def _settings_from_args(args):
    overrides = {name: getattr(args, name) for name in PARAMS}
    try:
        return build_settings(args.preset, args.config, overrides)
    except ValueError as e:
        raise SystemExit(f"error: {e}")


def _run(args):
    run_acquisition(args.input_folder, args.output, _settings_from_args(args))


def _entropy(args):
    settings = dict(DEFAULTS)
    if args.brightfield_suffix:
        settings["brightfield_suffix"] = args.brightfield_suffix
    df = entropy_table(args.input_folder, settings)
    if df.empty:
        raise SystemExit(f"no *{settings['brightfield_suffix']} images in {args.input_folder}")
    df.to_csv(args.output, index=False)
    print(df["entropy"].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).to_string())
    print(f"\nWritten to {args.output}. Inspect the extremes to choose --entropy-min / --entropy-max.")


def _demo(args):
    from .demo import make_demo_dataset

    folders = make_demo_dataset(args.output_folder, replicates=args.replicates)
    print(f"Demo images written to {', '.join(str(f) for f in folders)}")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="biosensor-quant", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="quantify single-cell biosensor ratios")
    run.add_argument("input_folder", help="folder with Cellpose masks and channel images")
    run.add_argument("-o", "--output", default="results", help="output folder [default: results]")
    run.add_argument("--preset", choices=sorted(PRESETS), help="channel preset")
    run.add_argument("--config", help="TOML file with settings (see config/example.toml)")
    _add_setting_options(run)
    run.set_defaults(func=_run)

    ent = sub.add_parser("entropy", help="tabulate brightfield entropy to choose a focus threshold")
    ent.add_argument("input_folder")
    ent.add_argument("-o", "--output", default="entropy.csv")
    ent.add_argument("--brightfield-suffix")
    ent.set_defaults(func=_entropy)

    demo = sub.add_parser("demo", help="write a synthetic demo dataset (no real images needed)")
    demo.add_argument("output_folder", nargs="?", default="demo")
    demo.add_argument("--replicates", type=int, default=3)
    demo.set_defaults(func=_demo)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
