#!/usr/bin/env python
# -*- coding: utf-8 -*-

import vtk


def main():
    namedColors = vtk.vtkNamedColors()

    uGrid = MakeQuadraticHexahedron()

    tessellate = vtk.vtkTessellatorFilter()
    tessellate.SetInputData(uGrid)
    tessellate.SetChordError(0.035)
    tessellate.Update()

    cellMap = dict()

    numTets = 0
    cell = vtk.vtkGenericCell()
    it = tessellate.GetOutput().NewCellIterator()
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        it.GetCell(cell)
        cellMap[cell.GetRepresentativeCell().GetClassName()] = numTets
        numTets += 1
        it.GoToNextCell()

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection(tessellate.GetOutputPort())
    mapper.ScalarVisibilityOff()

    # Create an actor for the grid
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetDiffuseColor(
        namedColors.GetColor3d("Tomato"))
    actor.GetProperty().SetEdgeColor(
        namedColors.GetColor3d("IvoryBlack"))
    actor.GetProperty().EdgeVisibilityOn()

    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetRadius(0.02)

    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetInputData(uGrid)
    glyph3D.SetSourceConnection(sphereSource.GetOutputPort())
    glyph3D.ScalingOff()
    glyph3D.Update()

    glyph3DMapper = vtk.vtkDataSetMapper()
    glyph3DMapper.SetInputConnection(glyph3D.GetOutputPort())
    glyph3DMapper.ScalarVisibilityOff()

    glyph3DActor = vtk.vtkActor()
    glyph3DActor.SetMapper(glyph3DMapper)
    glyph3DActor.GetProperty().SetColor(
        namedColors.GetColor3d("Banana"))

    textProperty = vtk.vtkTextProperty()
    textProperty.SetFontSize(24)

    ss = "# of Tetras: " + str(numTets)
    textMapper = vtk.vtkTextMapper()
    textMapper.SetInput(ss)
    textMapper.SetTextProperty(textProperty)

    textActor = vtk.vtkActor2D()
    textActor.SetMapper(textMapper)
    textActor.SetPosition(10, 400)

    # Visualize
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName("Quadratic Hexahedron Demo")
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(640, 512)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderWindow)

    widget = vtk.vtkSliderWidget()
    MakeWidget(widget, tessellate, textMapper, interactor)

    renderer.AddActor(actor)
    renderer.AddActor(glyph3DActor)
    renderer.AddViewProp(textActor)
    renderer.SetBackground(namedColors.GetColor3d("SlateGray"))

    renderWindow.Render()

    interactor.Start()


class SliderCallbackChordError():
    def __init__(self, tessellate, textMapper):
        self.tessellate = tessellate
        self.textMapper = textMapper

    def __call__(self, caller, ev):
        sliderWidget = caller
        value = sliderWidget.GetRepresentation().GetValue()
        self.tessellate.SetChordError(value)
        self.tessellate.SetMaximumNumberOfSubdivisions(4)
        self.tessellate.Update()

        cellMap = dict()

        numTets = 0
        cell = vtk.vtkGenericCell()
        it = self.tessellate.GetOutput().NewCellIterator()
        it.InitTraversal()
        while not it.IsDoneWithTraversal():
            it.GetCell(cell)
            cellMap[cell.GetRepresentativeCell().GetClassName()] = numTets
            numTets += 1
            it.GoToNextCell()
        ss = "# of Tetras: " + str(numTets)
        self.textMapper.SetInput(ss)


def MakeWidget(widget, tessellate, textMapper, interactor):
    # Setup a slider widget for each varying parameter
    tubeWidth = 0.008
    sliderLength = 0.008
    titleHeight = 0.04
    labelHeight = 0.04

    sliderRepChordError = vtk.vtkSliderRepresentation2D()

    sliderRepChordError.SetMinimumValue(0.0)
    sliderRepChordError.SetMaximumValue(0.07)
    sliderRepChordError.SetValue(tessellate.GetChordError())
    sliderRepChordError.SetTitleText("Chord error")

    sliderRepChordError.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    sliderRepChordError.GetPoint1Coordinate().SetValue(0.1, 0.1)
    sliderRepChordError.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    sliderRepChordError.GetPoint2Coordinate().SetValue(0.9, 0.1)

    sliderRepChordError.SetTubeWidth(tubeWidth)
    sliderRepChordError.SetSliderLength(sliderLength)
    sliderRepChordError.SetTitleHeight(titleHeight)
    sliderRepChordError.SetLabelHeight(labelHeight)

    widget.SetInteractor(interactor)
    widget.SetRepresentation(sliderRepChordError)
    widget.SetAnimationModeToAnimate()
    widget.EnabledOn()

    widget.AddObserver(vtk.vtkCommand.InteractionEvent, SliderCallbackChordError(tessellate, textMapper))


def MakeQuadraticHexahedron():
    aHexahedron = vtk.vtkQuadraticHexahedron()
    points = vtk.vtkPoints()

    #  pcoords = aHexahedron.GetParametricCoords()
    pcoords = (
        1.171002565324547,0.5200067577948213,1.0471337011831903,
        1.2385803952163401,0.5409110302503057,0.04963871457913571,
        0.5651829875035435,-0.1958436668685778,-0.011422300470955826,
        0.49760515761175017,-0.2167479393240621,0.9860726861330986,
        0.43481701249645655,1.1958436668685777,1.0114223004709557,
        0.5023948423882498,1.2167479393240621,0.013927313866901325,
        -0.17100256532454705,0.4799932422051787,-0.04713370118319016,
        -0.2385803952163401,0.45908896974969426,0.9503612854208643,
        1.2047914802704436,0.5304588940225635,0.5483862078811629,
        0.9018816913599419,0.17253368169086392,0.019108207054089943,
        0.5313940725576468,-0.20629580309632,0.48732519283107145,
        0.8343038614681486,0.1516294092353795,1.0166031936581446,
        0.4686059274423532,1.20629580309632,0.5126748071689285,
        0.1656961385318514,0.8483705907646205,-0.016603193658144444,
        -0.20479148027044358,0.4695411059774365,0.45161379211883707,
        0.0981183086400581,0.8274663183091361,0.98089179294591,
        0.8029097889105017,0.8579252123316996,1.029278000827073,
        0.870487618802295,0.878829484787184,0.031783014223018546,
        0.19709021108949826,0.1420747876683004,-0.029278000827072992,
        0.12951238119770497,0.12117051521281597,0.9682169857769815,
    )
    rng = vtk.vtkMinimalStandardRandomSequence()
    points.SetNumberOfPoints(aHexahedron.GetNumberOfPoints())
    rng.SetSeed(5070)  # for testing
    for i in range(0, aHexahedron.GetNumberOfPoints()):
        perturbation = [0.0] * 3
        for j in range(0, 3):
            rng.Next()
            #  perturbation[j] = rng.GetRangeValue(-0.1, 0.1)
            perturbation[j] = rng.GetRangeValue(0, 0)
        aHexahedron.GetPointIds().SetId(i, i)
        points.SetPoint(i, pcoords[3 * i] + perturbation[0],
                        pcoords[3 * i + 1] + perturbation[1],
                        pcoords[3 * i + 2] + perturbation[2])

    # Add the points and hexahedron to an unstructured grid
    uGrid = vtk.vtkUnstructuredGrid()
    uGrid.SetPoints(points)
    uGrid.InsertNextCell(aHexahedron.GetCellType(), aHexahedron.GetPointIds())

    return uGrid


if __name__ == '__main__':
    main()
