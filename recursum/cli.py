"""Command-line interface for RECURSUM."""

import argparse
import sys
from pathlib import Path


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="RECURSUM: Recurrence relation code generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  recursum generate              Generate all code
  recursum generate --output ./  Generate to specific directory
  recursum watch                 Watch for changes and auto-regenerate
  recursum list                  List available recurrences
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Generate command
    gen_parser = subparsers.add_parser('generate', help='Generate all code')
    gen_parser.add_argument(
        '--output', type=Path, default=Path.cwd(),
        help='Output directory (default: current directory)'
    )

    # Watch command
    watch_parser = subparsers.add_parser(
        'watch', help='Watch for changes and auto-regenerate'
    )
    watch_parser.add_argument(
        '--project-root', type=Path, default=Path.cwd(),
        help='Project root directory (default: current directory)'
    )

    # List command
    list_parser = subparsers.add_parser('list', help='List available recurrences')
    list_parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Show detailed information'
    )

    # Parse arguments
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 0

    # Execute command
    try:
        if args.command == 'generate':
            cmd_generate(args.output)
        elif args.command == 'watch':
            cmd_watch(args.project_root)
        elif args.command == 'list':
            cmd_list(args.verbose)
        else:
            parser.print_help()
            return 1
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        return 130
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        return 1

    return 0


def cmd_generate(output_dir: Path):
    """Execute generate command."""
    from .codegen.orchestrator import generate_all

    print(f"Generating code to: {output_dir}")
    generate_all(output_dir)
    print("\n✓ Generation complete!")


def cmd_watch(project_root: Path):
    """Execute watch command."""
    try:
        from .watch import watch
    except ImportError:
        print("Error: watchdog package not installed")
        print("Install with: pip install watchdog")
        sys.exit(1)

    print(f"Starting file watcher for: {project_root}")
    watch(project_root)


def cmd_list(verbose: bool = False):
    """Execute list command."""
    from .recurrences import get_all_recurrences

    print("Available recurrences:\n")

    recurrences = get_all_recurrences()

    total = 0
    for module, recs in recurrences.items():
        print(f"  {module}:")
        for rec in recs:
            indices_str = ', '.join(rec.indices)
            runtime_str = ', '.join(rec.runtime_vars) if rec.runtime_vars else "none"

            if verbose:
                print(f"    - {rec.name}")
                print(f"        Indices: {indices_str}")
                print(f"        Runtime vars: {runtime_str}")
                print(f"        Namespace: {rec.namespace}")
                print(f"        Max indices: {rec.max_indices}")
            else:
                print(f"    - {rec.name} (indices: {indices_str})")

        total += len(recs)
        print()

    print(f"Total: {total} recurrences across {len(recurrences)} modules")


if __name__ == "__main__":
    sys.exit(main())
