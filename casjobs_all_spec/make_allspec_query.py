#++++++++++++++++++++++++++
#
# TITLE: make_new_querey 
#
# PURPOSE: This function makes the
#          casjobs query used to get
#          the necessary data from casjobs.
#          
# INPUTS: stripe_82: if set to 1, the stripe 82 is searched
#         if set to 0, the DR7 is searched
#         default = 0
#
#         start_objid: the minimum objID used to extract 
#                      from the database. It is used 
#                      to extract from the database in chunks
#
# OUTPUTS: The query string used by casjobs
#
# PROGRAM CALLS: NONE
#
# BY: Alan Meert
#     Department of Physics and Astronomy
#     University of Pennsylvania
#
# FOR: Mariangela Bernardi
#      Department of Physics and Astronomy
#      University of Pennsylvania
#
# DATE: 9 FEB 2012
#
#-----------------------------------
import traceback

def start_query():
    """The first part of the query"""
    return """DECLARE @BRIGHT bigint SET @BRIGHT=dbo.fPhotoFlags('BRIGHT')
DECLARE @CHILD bigint SET @CHILD=dbo.fPhotoFLAGS('CHILD')
DECLARE @DEBLENDED_AS_PSF bigint 
    SET @DEBLENDED_AS_PSF=dbo.fPhotoFLAGS('DEBLENDED_AS_PSF')
DECLARE @EDGE bigint SET @EDGE=dbo.fPhotoFlags('EDGE')
DECLARE @SATURATED bigint SET @SATURATED=dbo.fPhotoFlags('SATURATED')
DECLARE @NODEBLEND bigint SET @NODEBLEND=dbo.fPhotoFlags('NODEBLEND')
DECLARE @bad_flags bigint SET
@bad_flags=(@SATURATED|@BRIGHT|@EDGE|@NODEBLEND|@CHILD|@DEBLENDED_AS_PSF)
SELECT 
"""

def mid_query():
    """The middle part of the query"""
    return """p.objid, (p.flags & @bad_flags)  as badflag,  
p.run,p.rerun,p.camCol,p.field,p.obj, s.specobjid, s.plate, s.mjd, s.fiberid,  
p.ra as ra_gal, p.dec as dec_gal, s.z, 
p.petroR50_g, p.petroR50_r, p.petroR50_i
p.petroMag_g, p.petroMag_r, p.petroMag_i,
p.extinction_g, p.extinction_r, p.extinction_i,
x.z as photoz, x.zErr as photoz_err,
x.kcorr_g,x.kcorr_r,x.kcorr_i
INTO
mydb.{tablename}
FROM
(photoobj as p LEFT OUTER JOIN SpecObj as s on p.objID = s.BestObjID) 
LEFT OUTER JOIN Photoz as x on x.objid =p.objid, chunk c,  field f, segment g
WHERE
g.segmentID = f.segmentID and
f.fieldID = p.fieldID and  c.chunkID=g.chunkID
""" 

def end_query():
    """The end of the query common to all query types"""
    return """ORDER BY p.objid\n"""

def catalog_query(job_info):
    """Takes a dictionary called job_info, with relevant target, chunksize, and starting point and returns the chunk of size chunksize containing galaxies with objid>starting point. Searches for galaxies with PetroMag_r-extinction_r between 14 and 17.77 that are part of the spectroscopic database"""

    out = start_query()
    out += "top {chunk}\n"
    out += mid_query()
    out += """and p.objid > {lastnum} and 
(p.petroMag_r - p.extinction_r) between 14.0 and 17.77 and p.type = 3 and s.specclass = 2\n""" 
    out += end_query()  

    try:
        out = out.format(**job_info)
    except KeyError:
        print """WARNING: Not all query info was supplied to the querey generator:
You must supply a limiting chunk size, a starting photoobjid, and table name"""
        traceback.print_exc()

    return out


def field_query(job_info):
    """returns a catalog of all objects in a given field"""

    out = start_query() + mid_query()
    out += """and p.run={run} and p.rerun = {rerun} and p.camcol = {camcol} 
and p.field = {field}\n""" 
    out += end_query()  

    try:
        out = out.format(**job_info)
    except KeyError:
        print """WARNING: Not all query info was supplied to the querey generator:
You must supply a run, rerun, camcol, field, and table name"""
        traceback.print_exc()

    return out

        
        


   
