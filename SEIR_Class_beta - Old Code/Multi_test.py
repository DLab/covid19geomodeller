#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import Multi_dist_ref_SEIR as MDSEIR

state='05'
comunas=['Valparaíso','Concón','Viña del Mar','Quilpué','Villa Alemana']
    

test=MDSEIR.refine_multi(state,comunas)