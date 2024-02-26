from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import pandas as pd


def data_from_single_tank(tank_id):
    # Filter the data for the tank_id
    tank_data = df[df["WasteSiteId"] == tank_id]
    return tank_data


app = Dash(__name__)

df = pd.read_csv("Tanks_Slurry_Inventory - all_tank_data.csv")

app.layout = html.Div(
    [
        html.H4("Tank analysis"),
        dcc.Graph(id="graph"),
        html.P("Tank ID:"),
        dcc.Dropdown(
            id="tankID",
            value="241-TX-101",
            options=df["WasteSiteId"].unique(),
        ),
        html.P("Values:"),
        dcc.Dropdown(
            id="values",
            options=["Mass (kg)", "Activity (Ci)", "WastePhase Mass (kg)"],
            value="Mass (kg)",
            clearable=False,
        ),
    ]
)


@app.callback(
    Output("graph", "figure"), Input("tankID", "value"), Input("values", "value")
)
def generate_chart(tankID, field):
    df_tank = data_from_single_tank(tankID)
    df_tank = df_tank.groupby("Analyte").sum()
    df_tank = df_tank.reset_index()
    print(df_tank)
    # only represent large values
    df_tank.loc[df_tank[field] < 1 / 100 * df_tank[field].sum(), "Analyte"] = "Others"

    fig = px.pie(df_tank, values=field, names="Analyte", hole=0.3)
    return fig


app.run_server(debug=True)
