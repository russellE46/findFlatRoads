import streamlit as st
import pandas as pd
import pydeck as pdk
import requests
import findRoads

tab = "‎ ‎ ‎ ‎ ‎ ‎ ‎ ‎ "
straightnessValues = {"Not Necessary": 10, "Kinda Straight": 0.2, "Straight": 0.08, "Very Straight": 0.01}
geo_api_key = "66048c06e66f8247638841qmoe4e996"

DATA_URL = "defaultRoadOutput.json"
if "roadData" not in st.session_state:
    df = pd.read_json(DATA_URL)

else:
    DATA_URL = "roadOutput.json"
    try:
        df = pd.read_json(st.session_state.roadData)
    except:
        df = df = pd.read_json(DATA_URL)

def get_coordinates(address, postal_code, api_key):
    # Base URL for the geocoding API
    base_url = "https://geocode.maps.co/search"

    # Parameters for the API request
    params = {"street": address, "postalcode": postal_code, "api_key": api_key}

    # Send GET request to the API
    response = requests.get(base_url, params=params)

    # Check if the request was successful
    if response.status_code == 200:
        # Parse the response JSON
        if len(response.json()) > 0:
            data = response.json()
            return (data[0]['lat'], data[0]['lon'])

    # If the request was unsuccessful, print an error message
    else:
        print("Error:", response.status_code)
    return (37.7009640, -121.9226350)

def hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

df['color'] = df['color'].apply(hex_to_rgb)

layer = pdk.Layer(
    type='PathLayer',
    data=df,
    pickable=True,
    get_color='color',
    width_scale=1,
    width_min_pixels=1,
    get_path='path',
    get_width=30
)

with st.sidebar:
    st.image("https://i.imgur.com/LRRnBkD.png", width = 200)
    #st.image("https://pbs.twimg.com/media/D5l75s8WwAIDDV8.jpg", width = 200)
    metricImp = st.selectbox(label="Units" , options=["Imperial","Metric"], index=0)

    if metricImp == "Metric":
        unitString = " (m)"
        radiusString = " (km)"
        defaultMinDistance = 400
        defaultMaxDistance = 460
        defaultMaxElevationChange = 4
        defaultMinRadius = 2.0
        defaultMaxRadius = 12.0
        defaultRadiusValue = 3.2
        rScale = 1050
        met = True

    else:
        unitString = " (ft)"
        radiusString = " (mi)"
        defaultMinDistance = 1320
        defaultMaxDistance = 1500
        defaultMaxElevationChange = 10
        defaultMinRadius = 1.0
        defaultMaxRadius = 8.0
        defaultRadiusValue = 1.8
        rScale = 1850
        met = False


    addressCoord = st.selectbox(label="Location Type" , options=["Address","Coordinates"], index=0)

    if addressCoord == "Address":
        addOld = ""
        postOld = ""
        if "address" not in st.session_state:
            st.session_state.address = st.text_input(label="Street Address", value="6568 Village Parkway")
            st.session_state.postal = st.text_input(label="Postal Code", value="94568")
        else:
            addOld = st.session_state.address
            postOld = st.session_state.postal
            st.session_state.address = st.text_input(label="Street Address", value=st.session_state.address)
            st.session_state.postal = st.text_input(label="Postal Code", value=st.session_state.postal)

        if addOld != st.session_state.address or postOld != st.session_state.postal:
            fetchedCoords = get_coordinates(st.session_state.address, st.session_state.postal, geo_api_key)
            st.session_state.origLat = float(fetchedCoords[0])
            st.session_state.origLon = float(fetchedCoords[1])

    else:
        if "origLat" not in st.session_state:
            st.session_state.origLat = st.number_input(label="Latitude:", min_value=-90.0, max_value=90.0, value=37.7009640, step=0.0000001, format="%.7f")
            st.session_state.origLon = st.number_input(label="Longitude:", min_value=-180.0, max_value=180.0, value=-121.9226350, step=0.0000001, format="%.7f") 
        else:
            st.session_state.origLat = st.number_input(label="Latitude:", min_value=-90.0, max_value=90.0, value=st.session_state.origLat, step=0.0000001, format="%.7f")
            st.session_state.origLon = st.number_input(label="Longitude:", min_value=-180.0, max_value=180.0, value=st.session_state.origLon, step=0.0000001, format="%.7f") 

    minLength = st.number_input(label="Minimum Distance" + unitString, min_value=100, max_value=2640, value=defaultMinDistance, step=20)
    maxLength = st.number_input(label="Maximum Distance" + unitString, min_value=150, max_value=3000, value=defaultMaxDistance, step=20)
    straightness = st.selectbox(label="Straightness", options=["Not Necessary","Kinda Straight", "Straight", "Very Straight"], index=2)
    if straightness == "Not Necessary":
        if metricImp == "Metric":
            defaultMaxRadius = 5.0
        else:
            defaultMaxRadius = 3.0
    elif straightness == "Kinda Straight":
        if metricImp == "Metric":
            defaultMaxRadius = 7.0
        else:
            defaultMaxRadius = 5.0
    elif straightness == "Straight":
        if metricImp == "Metric":
            defaultMaxRadius = 10.0
        else:
            defaultMaxRadius = 7.0
    maxElevChange = st.number_input(label="Maximum Elevation Change" + unitString, min_value=1, value=defaultMaxElevationChange, step=1)
    radius = st.slider(label="Radius" + radiusString, min_value=defaultMinRadius, max_value=defaultMaxRadius, value=defaultRadiusValue, step=0.1)
    st.caption("Note: maximum radius changes based on straightness requested")
if "pressed" not in st.session_state:
    st.title("Welcome!")
    st.markdown("Level is a tool to help you find flat and/or straight roads in your area.")
    st.markdown("Simply enter your location and desired road characteristics.")
    st.divider()
    st.markdown("Some notes:")
    st.markdown(tab + "• Residential roads will NOT be included in these searches.")
    st.markdown(tab + "• You can enter coordinates into Google maps, Apple Maps, etc.")
    st.markdown(tab + "• The more restrictive your search criteria is, the faster results will be returned.")
    st.markdown(tab + tab + "-- For example, requesting roads within a 2 mile radius will be much faster than requesting roads\n" + tab + tab + tab + "within a 10 mile radius")
    st.markdown(tab + "• In a similar vein, this tool relies on a free public API for elevation data. Because of this, we can only\n" + tab + tab + "process 100 coordinates per second. Please be mindful of this.")
    st.markdown(tab + "• Some road segments will be displayed as shorter than they actually are. This is a known bug.")
    st.markdown(tab + "• I am not responsible for what you do with the information provided by this website.")
    st.divider()
else:
    st.title("Results:")

origin =[
    {
        "name": str(st.session_state.origLat) + ", " + str(st.session_state.origLon),
        "color": "#ebd007",
        "coordinates": [
            st.session_state.origLon, st.session_state.origLat
        ]
    }
]
originDF = pd.DataFrame(origin)
radiusLayer = pdk.Layer(
    "ScatterplotLayer",
    data=originDF,
    pickable=False,
    opacity=0.008,
    stroked=True,
    filled=True,
    radius_scale=rScale,
    radius_min_pixels=1,
    radius_max_pixels=1000,
    line_width_min_pixels=3,
    get_position="coordinates",
    get_radius=radius,
    get_fill_color=[247, 114, 205],
    get_line_color=[255, 72, 192],
)

view_state = pdk.ViewState(
    latitude=st.session_state.origLat,     
    longitude=st.session_state.origLon,
    zoom=13,
    height = 850,
    width = "110%"
)


r = pdk.Deck(layers=[layer, radiusLayer], initial_view_state=view_state, tooltip={'text': '{name}'})
st.pydeck_chart(r)

if "roadDetails" in st.session_state:
    try:
        st.caption("Click on a column title to sort by that attribute")
        detailDf = pd.read_json(st.session_state.roadDetails)
        st.dataframe(data=detailDf, hide_index=True, use_container_width=True, height=500)
    except Exception as e:
        detailDf = []
else:
    try:
        st.caption("Click on a column title to sort by that attribute")
        detailDf = pd.read_json("defaultRoadDetailOutput.json")
        st.dataframe(data=detailDf, hide_index=True, use_container_width=True, height=500)
    except Exception as e:
        detailDf = []


st.sidebar.button(label="Find", on_click=findRoads.findRoads, args=(st.session_state.origLat,st.session_state.origLon,radius,minLength,
                                       maxLength,straightnessValues[straightness], maxElevChange, met))
