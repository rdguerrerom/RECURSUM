"""File watcher for automatic code regeneration."""

import time
import subprocess
from pathlib import Path

try:
    from watchdog.observers import Observer
    from watchdog.events import FileSystemEventHandler
    WATCHDOG_AVAILABLE = True
except ImportError:
    WATCHDOG_AVAILABLE = False
    Observer = None
    FileSystemEventHandler = object


class RecurrenceChangeHandler(FileSystemEventHandler):
    """Handle changes to recurrence definition files."""

    def __init__(self, project_root: Path):
        """
        Initialize change handler.

        Args:
            project_root: Root directory of the project
        """
        self.project_root = project_root
        self.last_rebuild = 0
        self.debounce_seconds = 2

    def on_modified(self, event):
        """Handle file modification events."""
        if event.is_directory:
            return

        # Only trigger on .py files in recurrences/ directory
        if not event.src_path.endswith('.py'):
            return
        if 'recurrences' not in event.src_path:
            return

        # Debounce rapid changes
        now = time.time()
        if now - self.last_rebuild < self.debounce_seconds:
            return

        self.last_rebuild = now
        print(f"\n📝 Detected change: {event.src_path}")
        self.rebuild()

    def rebuild(self):
        """Regenerate code and rebuild C++ extensions."""
        print("\n" + "=" * 70)
        print("🔄 Regenerating code...")
        print("=" * 70)

        try:
            # Regenerate code
            from .codegen.orchestrator import generate_all
            generate_all(self.project_root)
            print("\n✓ Code generation complete")

            # Rebuild C++ extensions
            print("\n" + "=" * 70)
            print("🔨 Rebuilding C++ extensions...")
            print("=" * 70)

            result = subprocess.run(
                ["pip", "install", "-e", ".", "--no-build-isolation"],
                cwd=self.project_root,
                capture_output=True,
                text=True
            )

            if result.returncode == 0:
                print("\n✓ Build complete! Extensions reloaded.")
                print("\nNote: Restart Python to use updated extensions.")
            else:
                print(f"\n✗ Build failed:")
                print(result.stderr)

        except Exception as e:
            print(f"\n✗ Error during rebuild: {e}")
            import traceback
            traceback.print_exc()


def watch(project_root: Path = None):
    """
    Start file watcher for auto-regeneration.

    Args:
        project_root: Root directory of the project (default: current directory)
    """
    if not WATCHDOG_AVAILABLE:
        print("Error: watchdog package not installed")
        print("Install with: pip install watchdog")
        return

    if project_root is None:
        project_root = Path.cwd()
    elif isinstance(project_root, str):
        project_root = Path(project_root)

    watch_dir = project_root / "recursum" / "recurrences"

    if not watch_dir.exists():
        print(f"Error: Directory not found: {watch_dir}")
        return

    print(f"👁️  Watching {watch_dir} for changes...")
    print("Press Ctrl+C to stop.\n")

    event_handler = RecurrenceChangeHandler(project_root)
    observer = Observer()
    observer.schedule(
        event_handler,
        str(watch_dir),
        recursive=True
    )
    observer.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
        print("\n\n👋 Watcher stopped.")

    observer.join()


if __name__ == "__main__":
    watch()
