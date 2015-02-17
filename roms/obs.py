#!/usr/bin/env python
"""
  obs.py
  
  State Estimation and Analysis for PYthon

  Module to handle the observation structure within ROMS

  Written by Brian Powell on 08/05/14
  Copyright (c)2014 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np
from scipy import ndimage
import os
import re

# Define the observation provenances used within my applications
obs_provenance = {
  100:"GLIDER",
  114:"GLIDER_SG114",
  139:"GLIDER_SG139",
  146:"GLIDER_SG146",
  147:"GLIDER_SG147",
  148:"GLIDER_SG148",
  150:"GLIDER_SG500",
  200:"CTD",
  210:"CTD_HOT",
  220:"CTD_ARGO",
  300:"SST",
  301:"SST_OSTIA",
  317:"SST_AVHRR_17",
  318:"SST_AVHRR_18",
  330:"SST_MODIS_AQUA",
  331:"SST_MODIS_TERRA",
  400:"SSH",
  411:"SSH_AVISO_ENVISAT",
  412:"SSH_AVISO_JASON1",
  413:"SSH_AVISO_JASON2",
  414:"SSH_AVISO_GFO",
  415:"SSH_HYCOM",
  500:"DRIFTERS",
  600:"RADAR",
  610:"RADAR_KOKOHEAD",
  620:"RADAR_KAKAAKO",
  630:"RADAR_KALAELOA",
  700:"ADCP",
  711:"ADCP_KILOMOANA_os38nb",
  712:"ADCP_KILOMOANA_os38bb",
  713:"ADCP_KILOMOANA_os75nb",
  714:"ADCP_KILOMOANA_os75bb",
  715:"ADCP_KILOMOANA_wh300",
  721:"ADCP_WECOMA_os38nb",
  722:"ADCP_WECOMA_os38bb",
  723:"ADCP_WECOMA_os75nb",
  724:"ADCP_WECOMA_os75bb",
  725:"ADCP_WECOMA_wh300",
  731:"ADCP_KOK_os38nb",
  732:"ADCP_KOK_os38bb",
  733:"ADCP_KOK_os75nb",
  734:"ADCP_KOK_os75bb",
  735:"ADCP_KOK_wh300",
  741:"ADCP_HIIALAKAI_os38nb",
  742:"ADCP_HIIALAKAI_os38bb",
  743:"ADCP_HIIALAKAI_os75nb",
  744:"ADCP_HIIALAKAI_os75bb",
  745:"ADCP_HIIALAKAI_wh300",
  751:"ADCP_THOMPSON_os38nb",
  752:"ADCP_THOMPSON_os38bb",
  753:"ADCP_THOMPSON_os75nb",
  754:"ADCP_THOMPSON_os75bb",
  755:"ADCP_THOMPSON_wh300"
}