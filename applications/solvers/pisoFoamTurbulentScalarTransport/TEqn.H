// Solve the turbulent scalar transport equation
if(solveScalar)
{
    

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(DT, T)
      - fvm::laplacian(
                        turbulence->nut()/Sct,
                        T
                       )
      ==
        fvOptions(T)
    );

    TEqn.relax();
    fvOptions.constrain(TEqn);
    TEqn.solve();
    fvOptions.correct(T);

}
