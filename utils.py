"""
    Utility functions to plot the model history, fetch data on Google Earth Engine and more
"""
import pandas as pd
import plotly_express as px
import ee
from datetime import datetime as dt, timedelta
def plot_model_history(history_df: pd.DataFrame, variables_to_plot: list= [ 'soc', 'doc', 'mic', 'enz', 'co2']):
    """ plots selected variables over time in a line graph"""
    
    fig = px.line(history_df, x='index', y=variables_to_plot, facet_row='variable')

    fig.update_yaxes(matches=None, title="")
    fig.update_xaxes(title="")
    fig.update_layout(
            title = f"System evolution over time. Total runs: {len(history_df)}",
            geo=dict(bgcolor= 'rgba(255,255,255,255)'),
            paper_bgcolor = 'rgba(255,255,255,255)',
            plot_bgcolor = 'rgba(255,255,255,255)',
            font_color = 'rgba(0,0,0,255)'
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.add_annotation(x=0.5,y= -0.3,
                    text="time", textangle=-0,
                        xref="paper", yref="paper")
    # Show the figure
    return fig

def authenticate_ee(
        service_account = 'digitsoiltemp@deep-bolt-386711.iam.gserviceaccount.com',
        path_to_keys = 'deep-bolt-386711-f24058c11291.json'
):
    credentials = ee.ServiceAccountCredentials(service_account, path_to_keys )
    ee.Initialize(credentials)

def sample_image(coords:tuple, image_id:str, band_name:str):
    # Authenticate to GEE
    authenticate_ee()

    # Define the point of interest
    point = ee.Geometry.Point(*coords)

    # Load the image you want to sample
    image = ee.Image(image_id)

    # Sample the image at the point of interest
    value = image.sample(point, 30).first().get(band_name).getInfo()

    return value

# Load the image collection and get the most recent image
def sample_collection(coords, collection_id, band_name:str, start_date:str, end_date:str):
    # Authenticate to GEE
    authenticate_ee()

    # Define the point of interest
    point = ee.Geometry.Point(*coords)

    # Load the image collection and filter it by date
    collection = ee.ImageCollection(collection_id) \
        .filterDate(start_date, end_date) \
        .sort('system:time_start')

    # Get the most recent image in the collection
    image = ee.Image(collection.first())

    # Sample the image at the point of interest
    value = image.sample(point, 30).first().get(band_name).getInfo()

    return value

# get SOC values
def get_current_soc_value(selected_point: tuple)->str:
    soc_value = sample_image(selected_point, "OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02", "b10")
    #convert from 5*g/kg to mg/cm³. Average density of soil: 2.65 g/cm³
    soc_value_translated = soc_value /5 * 1000 /2.65
    return soc_value_translated

def get_soil_temp(selected_point):
    now = dt.now()

    # Calculate the date 20 days ago
    start = now - timedelta(days=20)

    temp_value = sample_collection(selected_point, 'MODIS/061/MOD11A2', 'LST_Day_1km', start.strftime('%Y-%m-%d'),now.strftime('%Y-%m-%d') )
    # convert from 8 bit kelvin to °C
    temp_celsius = temp_value * 0.02 - 273.15
    return temp_celsius

def get_temp_projection(selected_point):
    now = dt.now()

    # Calculate the date 5 days ago
    future = now + timedelta(days=60)

    temp_value = sample_collection(selected_point, "NASA/GDDP-CMIP6", 'tasmax', now.strftime('%Y-%m-%d'),future.strftime('%Y-%m-%d') )
    # convert from 8 bit kelvin to °C
    temp_celsius = temp_value - 273.15
    return temp_celsius
