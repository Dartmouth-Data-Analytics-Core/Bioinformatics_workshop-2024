# GitHub Pages Setup Guide

This repository is configured to automatically deploy to GitHub Pages using GitHub Actions.

## Automatic Deployment

The site will automatically build and deploy when you push to the `main`, `master`, or `bookdown` branch.

### Setup Steps:

1. **Enable GitHub Pages in your repository:**
   - Go to your repository on GitHub
   - Click on **Settings** → **Pages**
   - Under "Source", select **GitHub Actions** (not "Deploy from a branch")
   - Save the settings

2. **Push your changes:**
   - The GitHub Actions workflow (`.github/workflows/deploy-pages.yml`) will automatically:
     - Build the site using the `build_site.py` script
     - Deploy it to GitHub Pages
   - Your site will be available at: `https://[your-username].github.io/Bioinformatics_workshop-2024/`

## Manual Deployment (Alternative)

If you prefer to deploy manually:

1. Build the site locally:
   ```bash
   python3 build_site.py
   ```

2. Push the `docs/` directory to the `gh-pages` branch:
   ```bash
   git subtree push --prefix docs origin gh-pages
   ```

   Or if the `gh-pages` branch doesn't exist:
   ```bash
   git checkout --orphan gh-pages
   git rm -rf .
   cp -r docs/* .
   git add .
   git commit -m "Deploy site"
   git push origin gh-pages
   ```

3. Configure GitHub Pages to use the `gh-pages` branch as the source.

## Custom Domain (Optional)

If you want to use a custom domain:

1. Create a `CNAME` file in the `docs/` directory with your domain name
2. Update your DNS settings to point to GitHub Pages
3. The workflow will automatically include the CNAME file in deployments

## Troubleshooting

- **Images not showing:** Make sure all image paths in markdown files use `../figures/` or `figures/` (the build script normalizes these)
- **Build fails:** Check that Python 3.9+ is available and the `markdown` package is installed
- **Pages not updating:** Wait a few minutes after pushing, GitHub Pages can take 1-2 minutes to update

