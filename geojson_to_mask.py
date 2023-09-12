import json
import sys
import os
import cv2
import numpy as np
from PIL import Image
from skimage.io import imsave

from io import BytesIO
import requests

sys.path.append('../image_computing_base_bb')
from utils import APIutils, parseargs

def hex_to_rgb(hexa):
    return tuple(int(hexa[i:i+2], 16)  for i in (0, 2, 4))

def apply_correction(pts,rot,jp2_wid,jp2_hei):
    theta = 0
    if rot==270:
        theta = -90*np.pi/180
    elif rot==90:
        theta = 90*np.pi/180
    elif rot==180:
        theta = np.pi
#     elif rot == 0:
#         theta = 0

    c,s = np.cos(theta), np.sin(theta)

    # x' = x c - y s
    # y' = x s + y c

    org = np.array([jp2_wid/2,-jp2_hei/2])

    pts_rot = []
    for pt in pts:
        x,y = pt-org
        x_ = x*c-y*s+org[0]
        y_ = x*s+y*c+org[1]
        pts_rot.append([x_,y_])

    pts_rot = np.array(pts_rot)
    return pts_rot 

def fill_canvasFinal(data,rotation,width,height, ratio=32):

    mask_by_id = {} # id:mask
    name_by_id = {} # id:name
    color_by_id = {} # id:(r,g,b)

    blank_img = np.zeros((height//ratio, width//ratio),np.uint8)
    
    for feature in data['features']:
        if 'Polygon' not in feature['geometry']['type']:
            continue
        parts = feature['geometry']['coordinates']

        color = feature['properties']['data']['color_hex_triplet']
        acro = feature['properties']['data']['acronym']
        name = feature['properties']['data']['name'].replace(' ',"_")
        ontoid = feature['properties']['data']['id']
        color = hex_to_rgb(color)

        if ontoid not in mask_by_id:
            mask_by_id[ontoid]=blank_img.copy()
            name_by_id[ontoid]=f'{name}({acro})'
            color_by_id[ontoid]=color
        
        for part in parts:
            part = apply_correction(part, rotation, width, height)
            outer_np = np.squeeze(np.array(part))
            outer_np = np.vstack((outer_np))
            outer_np/=(ratio)
            outer_np = np.abs(outer_np).astype(int)
            # image_580 = np.ascontiguousarray(blank_img, dtype=np.uint8)
            cv2.fillPoly(mask_by_id[ontoid], [outer_np],255 )

    return mask_by_id, name_by_id, color_by_id

def get_geojson(bsid, secno):
    gjurllistfile = './geojson_urls_%d.json' % bsid
    geojson = None
    if os.path.exists(gjurllistfile):
        gjdata = json.load(open(gjurllistfile,'rt'))
        for elt in gjdata:
            if elt['id']==secno:
                geojson = requests.get(elt['url']+'?format=json',auth=('admin','admin')).json()
                break
    return json.loads(geojson['atlasgeojson'])

if __name__=="__main__":

    args = parseargs(
        {
            1:('biosampleid',int,147),
            2:('section',int,301),
            3:('ontologyid',int,6),
            4:('outdir',str,'./structuremasks')
        }
    )

    print(args)

    obj = APIutils()
    ssid = obj.get_ssid(args.biosampleid)
    imglist = obj.get_imglist_api(args.biosampleid)

    imginfo = imglist['NISL'][args.section]
    
    # ontologydata = obj.get_atlastree(args.ontologyid)

    # Subpial Granular Layer (SG) : 11580
    #  Marginal Zone(MZ) : 10508
    #  Cortical Plate (CP) : 10515
    #  Subplate (SP) : 10522
    #  Intermediate Zone IZ) : 10529
    #  Subventricular Zone (SVZ) : 10536
    #  Ventricular zone (VZ) : 10542

    bn = os.path.basename(imginfo['jp2Path']).replace('compressed.jp2','.json')
    
    geojson = get_geojson(args.biosampleid, args.section)
    assert geojson is not None
    # geojson = json.load(open(gj_dir+'/'+bn))

    ratio = 32 # = 16/0.5 ; to bring to 16mpp from 0.5mpp
    msks, names, colors = fill_canvasFinal(geojson, imginfo['rigidrotation'], imginfo['width'], imginfo['height'], ratio=32)

    
    resp = requests.get(f'http://apollo2.humanbrain.in:9081/fcgi-bin/iipsrv.fcgi?FIF={imginfo["jp2Path"]}&WID={imginfo["width"]//ratio}&CVT=jpeg')
    nisslimg = Image.open(BytesIO(resp.content))

    outloc = f'{args.outdir}/bs-{args.biosampleid}'
    os.makedirs(outloc,exist_ok=True)

    imsave(f'{outloc}/{args.section}_nissl.png',np.array(nisslimg))
    os.makedirs(f'{outloc}/{args.section}',exist_ok=True)
    json.dump(names,open(f'{outloc}/{args.section}/names.json','wt'))
    json.dump(colors,open(f'{outloc}/{args.section}/colors.json','wt'))

    for ontoid in msks:
        imsave(f'{outloc}/{args.section}/{ontoid}.png',msks[ontoid])

