""" cnvpytor.report

Class Report: pdf template for reports
"""
from __future__ import absolute_import, print_function, division

from fpdf import FPDF, XPos, YPos
import logging
import pkg_resources

_logger = logging.getLogger("cnvpytor.report")



class Report(FPDF):
    def __init__(self):
        FPDF.__init__(self)
        self.sample = "<SAMPLE_NAME>"
        self.set_title("CNVpytor report")
        self.set_author("CNVpytor")
        self.date = "1/1/2001"

    def load_resource(self, reason, filename):
        if reason == "image":
            if filename.startswith("http://") or filename.startswith("https://"):
                f = BytesIO(urlopen(filename).read())
            elif filename.startswith("data"):
                f = filename.split('base64,')[1]
                f = base64.b64decode(f)
                f = io.BytesIO(f)
            else:
                f = open(filename, "rb")
            return f
        else:
            self.error("Unknown resource loading reason \"%s\"" % reason)

    def header(self):
        # Colors of frame, background and text
        self.set_draw_color(0, 150, 0)
        self.set_fill_color(230, 230, 230)
        self.set_text_color(0, 0, 0)
        self.set_line_width(0.2)
        self.rect(x=9, y=2, w=self.epw, h=22, style="FD")
        self.set_font('Arial', 'B', 14)
        self.set_y(4)
        self.set_x(100)
        self.cell(self.epw-105, 6, "SINGLE SAMPLE REPORT", new_x=XPos.LEFT, new_y=YPos.NEXT)
        self.cell(self.epw-105, 6, "SAMPLE: {:}".format(self.sample), new_x=XPos.LEFT, new_y=YPos.NEXT)
        self.cell(self.epw-105, 6, "DATE: {:}".format(self.date), new_x=XPos.LEFT, new_y=YPos.NEXT)

        self.image(pkg_resources.resource_filename('cnvpytor', 'imgs')+"/cnvpytor_640.png",x=11,y=2,w=80)

        # Line break
        self.ln(5)

    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Text color in gray
        self.set_text_color(128)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'C')

    def add_title(self, title):
        # Arial 12
        self.set_font('Arial', '', 12)
        # Background color
        self.set_fill_color(200, 220, 255)
        # Title
        self.cell(0, 6, title, 0, 1, 'L', 1)
        # Line break
        self.ln(4)

    def add_paragraph(self, txt):
        # Times 12
        self.set_font('Times', '', 12)
        # Output justified text
        self.multi_cell(0, 5, txt)
        # Line break
        self.ln()

    def add_plot(self, img):
        self.image(img, w=self.epw)
        self.ln()
