#!/usr/bin/env python3

import pandas as pd
import numpy as np
from PIL import Image
from PIL import ImageDraw, ImageFont
import argparse as arp



prs = arp.ArgumentParser()

prs.add_argument('-im','--images',
                  type = str,
                  nargs = '+',
                  required = True,
                  help = '',
                  )

prs.add_argument('-bw','--border_width',
                  type = int,
                  required = False,
                  default = 0,
                  help = '',
                 )

prs.add_argument('-ah','--add_header',
                  default = False,
                  action = 'store_true',
                  help = '',
                 )

prs.add_argument('-hc','--header_color',
                 default = 'gray',
                 type = str,
                 help = '',
                )

prs.add_argument('-t','--title',
                 default = '',
                 type = str,
                 help = '',
                )

prs.add_argument('-f','--font_path',
                 default = None,
                 type = str,
                 help = '',
                )

prs.add_argument('-o','--output_name',
                 default = '/tmp/joined.png',
                 type = str,
                 help = '',
                )

prs.add_argument('-bc','--border_color',
                 default = 'white',
                 type = str,
                 help = '',
                )

args = prs.parse_args()



pths = (args.images if isinstance(args.images,list) else [args.images])
add_header = args.add_header

if args.font_path is None:
    font_pth = "/usr/share/fonts/cantarell/Cantarell-Regular.otf"
else:
    font_pth = args.font_path

images = list(map(Image.open,pths))
widths, heights = zip(*(i.size for i in images))
border_width = args.border_width
bc_color = args.border_color

total_width = sum(widths) + border_width*(len(widths) + 1)
if add_header:
    txt = args.title
    header_extra = int(images[0].size[1] / 5)
    header = Image.new('RGB',
                       size = (total_width,header_extra),
                       color = args.header_color,
                      )

    img_fraction = 0.50
    fontsize = 1
    font = ImageFont.truetype(font_pth, fontsize)
    while font.getsize(txt)[0] < img_fraction*total_width:
        fontsize += 1
        font = ImageFont.truetype(font_pth, fontsize)
else:
    header_extra = 0

max_height = max(heights)
print(max_height,total_width)

new_im = Image.new('RGB',
                   size = (total_width, max_height + header_extra),
                   color = bc_color
                  )

x_offset = 0
x_offset += border_width

for im in images:
    new_im.paste(im, (x_offset,header_extra))
    x_offset += im.size[0] + border_width
if add_header:
    new_im.paste(header,(0,0))
    draw = ImageDraw.Draw(new_im)
    w,h = draw.textsize(txt, font = font)
    draw.text(((total_width -w)/2, h/2),
              txt,
              (0, 0, 0),
              font = font,
              )

new_im.save(args.output_name)
