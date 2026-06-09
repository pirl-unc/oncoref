# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""``cancerdata`` command-line interface.

Reference lookups over the bundled tables (cancer-type / TMB / burden) plus the
data cache surface. The heavy per-cohort expression bundle and its
``fetch``/``status``/``prune`` subcommands are added in a later milestone; the
cache-dir resolution here is already the path that bundle will populate.
"""

from __future__ import annotations

import argparse
import json
import sys

from . import cancer_types, incidence, tmb
from .cache import bundle_cache_dir
from .version import __version__


def _cmd_version(args: argparse.Namespace) -> int:
    print(f"cancerdata v{__version__}")
    return 0


def _cmd_cache_dir(args: argparse.Namespace) -> int:
    print(bundle_cache_dir())
    return 0


def _cmd_cancer_type(args: argparse.Namespace) -> int:
    try:
        info = cancer_types.cancer_type_info(args.query)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    if info is None:
        print("Error: no cancer type given", file=sys.stderr)
        return 1
    print(json.dumps(info, indent=2, default=str))
    return 0


def _cmd_tmb(args: argparse.Namespace) -> int:
    if args.code is None:
        for code, value in sorted(tmb.cancer_tmb().items()):
            print(f"{code}\t{value:g}")
        return 0
    try:
        value = tmb.cancer_tmb(args.code, inherit=not args.no_inherit)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    if value is None:
        print(f"No TMB value for {args.code!r}", file=sys.stderr)
        return 1
    print(f"{value:g}")
    return 0


def _cmd_burden(args: argparse.Namespace) -> int:
    try:
        if args.category is None:
            for cat, pct in sorted(incidence.cancer_burden(metric=args.metric).items()):
                print(f"{cat}\t{pct:g}")
            return 0
        value = incidence.cancer_burden(args.category, metric=args.metric)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    if value is None:
        print(f"No burden value for category {args.category!r}", file=sys.stderr)
        return 1
    print(f"{value:g}")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cancerdata",
        description="Curated cancer reference data: ontology, TMB, incidence/mortality.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("version", help="Print the installed cancerdata version").set_defaults(
        func=_cmd_version
    )

    sub.add_parser(
        "cache-dir", help="Print the on-disk cache dir for the downloadable data bundle"
    ).set_defaults(func=_cmd_cache_dir)

    p_ct = sub.add_parser(
        "cancer-type", help="Resolve a cancer type/alias/name and print its registry info"
    )
    p_ct.add_argument("query", help="Cancer code, alias, or display name (e.g. PRAD, prostate)")
    p_ct.set_defaults(func=_cmd_cancer_type)

    p_tmb = sub.add_parser("tmb", help="Median TMB (mut/Mb) for a code, or the full map")
    p_tmb.add_argument("code", nargs="?", default=None, help="Cancer code/alias (omit for all)")
    p_tmb.add_argument("--no-inherit", action="store_true", help="Do not inherit an ancestor's TMB")
    p_tmb.set_defaults(func=_cmd_tmb)

    p_burden = sub.add_parser(
        "burden", help="Incidence/mortality share for a burden category, or the full map"
    )
    p_burden.add_argument(
        "category", nargs="?", default=None, help="Burden category (omit for all)"
    )
    p_burden.add_argument(
        "--metric",
        default="us_incidence_pct",
        choices=list(incidence._BURDEN_METRICS),
        help="Which share to report (default: us_incidence_pct)",
    )
    p_burden.set_defaults(func=_cmd_burden)

    return parser


def main(argv=None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
