"""Entry point: launch the interactive airfoil flow explorer.

    uv run python src/main.py        # launch the interactive window
    uv run python src/validate.py    # run the physics validation suite
"""

from aerosim.app import main

if __name__ == "__main__":
    main()
