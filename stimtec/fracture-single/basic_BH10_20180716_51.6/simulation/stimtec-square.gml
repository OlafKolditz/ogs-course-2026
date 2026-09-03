<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>stimtec</name>
    <points>
        <point id="0" x="-50.0" y="-50.0" z="0" name="inflow"/>
        <point id="1" x="50.0" y="-50.0" z="0"/>
        <point id="2" x="50.0" y="50.0" z="0" name="outflow"/>
        <point id="3" x="-50.0"  y="50.0" z="0"/>
        <point id="4" x="0"  y="0" z="0" name="well"/>
    </points>
    <polylines>
        <polyline id="0" name="left">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>		
        <polyline id="1" name="right">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>		
        <polyline id="2" name="boundary">
            <pnt>0</pnt>
            <pnt>1</pnt>
            <pnt>2</pnt>
            <pnt>3</pnt>
        </polyline>	
    </polylines>
</OpenGeoSysGLI>