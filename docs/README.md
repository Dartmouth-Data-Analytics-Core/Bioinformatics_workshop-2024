# Static Website

This directory contains the generated static website for the Fundamentals of Bioinformatics Workshop.

## Viewing the Site

Simply open `index.html` in your web browser. You can:

1. Double-click the `index.html` file in your file manager
2. Or open it from the command line:
   ```bash
   open index.html  # macOS
   # or
   xdg-open index.html  # Linux
   # or
   start index.html  # Windows
   ```

## Rebuilding the Site

If you make changes to the markdown files, rebuild the site by running:

```bash
python3 build_site.py
```

This will regenerate the `index.html` file with your latest changes.

## Deployment

The `docs/` directory contains all the files needed for the static website. You can:

- Upload the entire `docs/` directory to any web hosting service
- Deploy to GitHub Pages (see `GITHUB_PAGES_SETUP.md` in the root directory)
- Use any static site hosting service (Netlify, Vercel, etc.)

### GitHub Pages

This site is configured for automatic deployment to GitHub Pages via GitHub Actions. See `GITHUB_PAGES_SETUP.md` in the repository root for detailed setup instructions.

The `.nojekyll` file in this directory tells GitHub Pages not to process the site with Jekyll, which is required for proper rendering.

