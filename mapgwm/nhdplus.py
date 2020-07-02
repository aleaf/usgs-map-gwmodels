import os
import gisutils


def get_flowline_routing(plusflow_file, dest_routing_file):
    if not os.path.exists(dest_routing_file):
        df = gisutils.shp2df(plusflow_file)
        routing = df[['FROMCOMID', 'TOCOMID']]
        routing.to_csv(dest_routing_file, index=False)