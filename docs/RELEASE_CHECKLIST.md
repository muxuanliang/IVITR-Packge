# Pre-release checklist

Before making this directory public:

1. Replace `TODO@example.com` in `code/IVITR/DESCRIPTION` with the package
   maintainer's public email address.
2. Confirm the application-data access statement in
   `docs/DATA_AVAILABILITY.md` with the PCPT data owner. If the analytic data
   can be publicly deposited, add its stable DOI and update that document.
3. Create a public GitHub repository from this directory, tag the exact
   manuscript version, and archive that tag with Zenodo (or another DOI-minting
   repository).
4. Insert the GitHub URL and version-specific DOI in the manuscript's data and
   code availability statement.
5. Run the commands in `docs/REPRODUCIBILITY.md` in a clean R installation.
