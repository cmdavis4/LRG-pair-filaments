from astro_image_processing.user_settings import casjobs_info

gal_cat = {'filename':'spectro_sample_raw.cat',
           'data_dir':'./data/',
           'out_file':'spectro_sample.cat'
           }

casjobs_info.update({ 'cas_jar_path':'/home/alan/git_projects/astro_image_processing/casjobs_query/casjobs.jar',
                 'jobname':'allspec',
                 'search_target':'DR10'})
