""" cnvpytor.annotator

Annotator class

"""

from __future__ import absolute_import, print_function, division
import requests
from .genome import *
from .utils import decode_region

_logger = logging.getLogger("cnvpytor.annotator")


class Annotator:

    def __init__(self, reference_genome, biotypes=["protein_coding"]):
        """
        Use Ensembl API rest to get annotation information for regions.

        Parameters
        ----------
        reference_genome : str
            Reference genome used. Expect reference genome configuration with field "ensembl_api_region"
        biotypes : list of str
            List of biotypes to annotate
        """
        self.reference_genome = reference_genome
        self.url_template = reference_genome["ensembl_api_region"]
        self.biotypes = biotypes

    def get_info(self, region):
        """
        Returns annotation
        Parameters
        ----------
        region : str
            Region in genome

        Returns
        -------
        ret : string
            Annotation string
        """
        regs = decode_region(region)
        if len(regs) == 0:
            _logger.debug("get_info called without region")
            return None
        (c, (start, end)) = regs[0]
        response = requests.get(self.url_template.format(region=region))
        ret = []
        if "error" in response.json():
            return "No information"
        for i in response.json():
            if ("biotype" in i) and i["biotype"] in self.biotypes:
                if i["start"] > start and i["end"] < end:
                    position = "inside"
                elif i["start"] < start and i["end"] > end:
                    position = "cover"
                elif i["start"] < start:
                    position = "intersect left"
                else:
                    position = "intersect right"
                if "external_name" in i:
                    ret.append("%s (%s %s)" % (i["external_name"], i["id"], position))
                else:
                    ret.append("%s (%s %s)" % (i["id"], i["id"], position))
        return ", ".join(ret)
