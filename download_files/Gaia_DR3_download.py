
import os
import gzip
import requests
from bs4 import BeautifulSoup
import pandas as pd
import time as t


gaia_url = 'http://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/'


def main(mag_max=19, dest_folder='datafiles'):
    """
    Download Gaia DR3 source data files.

    Description of files:
    https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html
    """

    gaia_files_file = "Gaia_DR3_files.txt"

    if os.path.exists(gaia_files_file):
        # Read file
        with open(gaia_files_file, 'r') as file:
            gaia_files = file.readlines()
        gaia_files = [_.strip() for _ in gaia_files]
        print("File with Gaia DR3 list of files read")
    else:
        # Generate file
        ext = 'gz'
        gaia_files = get_url_paths(gaia_url, ext)

        with open(gaia_files_file, 'w') as file:
            for line in gaia_files:
                file.write(line + '\n')
        print("File with Gaia DR3 list of files generated")

    txt_f = 'frame_ranges.txt'
    if not os.path.exists(txt_f):
        with open(txt_f, 'w') as file:
            file.write('filename,ra_min,ra_max,dec_min,dec_max\n')

    req_cols = (
        'source_id', 'ra', 'dec', 'l', 'b',
        'parallax', 'parallax_error',
        'pmra', 'pmra_error', 'pmdec', 'pmdec_error',
        'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_mag',
        'phot_bp_mean_flux', 'phot_bp_mean_flux_error',
        'phot_rp_mean_flux', 'phot_rp_mean_flux_error',
        'radial_velocity', 'radial_velocity_error')

    for file_url in gaia_files:
        start = t.time()
        filename, file_path = download(file_url, dest_folder=dest_folder)

        try:
            with gzip.open(file_path) as f:
                data = pd.read_csv(
                    f, usecols=req_cols, index_col=False, comment='#')

                msk = data['phot_g_mean_mag'] < mag_max
                data = data[msk]

                # Remove to save space
                data = data.drop(columns=['phot_g_mean_mag'])

                ra_r = "{:.3f},{:.3f}".format(
                    data['ra'].min(), data['ra'].max())
                de_r = "{:.3f},{:.3f}".format(
                    data['dec'].min(), data['dec'].max())
                txt = filename.replace(
                    '.csv.gz', '_p.csv.gz') + ',' + ra_r + ',' + de_r

                with open(txt_f, 'a') as file:
                    file.write(txt + '\n')

                data.to_csv(
                    file_path.replace('.csv.gz', '_p.csv.gz'), index=False,
                    compression='gzip')

            os.remove(file_path)
        except Exception as err:
            print(err)

        print(t.time() - start)


def get_url_paths(url, ext='', params={}):
    """
    https://stackoverflow.com/a/49567481/1391441
    """
    response = requests.get(url, params=params)
    if response.ok:
        response_text = response.text
    else:
        return response.raise_for_status()
    soup = BeautifulSoup(response_text, 'html.parser')
    parent = [url + node.get('href') for node in soup.find_all('a')
              if node.get('href').endswith(ext)]
    return parent


def download(url: str, dest_folder: str):
    """
    https://stackoverflow.com/a/56951135/1391441
    """
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)  # create folder if it does not exist

    # be careful with file names
    filename = url.split('/')[-1].replace(" ", "_")
    file_path = os.path.join(dest_folder, filename)

    r = requests.get(url, stream=True)
    if r.ok:
        print("saving to", os.path.abspath(file_path))
        with open(file_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024 * 8):
                if chunk:
                    f.write(chunk)
                    f.flush()
                    os.fsync(f.fileno())
    else:  # HTTP status code 4XX/5XX
        print("Download failed: status code {}\n{}".format(
            r.status_code, r.text))

    return filename, file_path


if __name__ == '__main__':
    main()
