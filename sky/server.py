"""Local HTTP access to the same API used by the CLI."""

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from .api import RunContext, archive_call, cone, orbital_decay, survey_names
from .archives import SURVEYS, gaia_epoch_photometry, jwst_products


class RunRequest(BaseModel):
    surveys: str = "ztf,gaia,roman"
    top: int = Field(default=10, ge=1, le=100)
    limit_per_partition: int | None = Field(default=None, ge=1)


def create_app(partitions=(), *, workers=1):
    app = FastAPI(title="Sky survey API", version="0.1.0")

    @app.get("/v1/surveys")
    def surveys():
        return SURVEYS

    @app.get("/v1/cone")
    def cone_search(ra: float, dec: float, radius_arcsec: float = 60, limit: int = 20,
                    surveys: str = "gaia,roman", roman_program: str = "gbtds"):
        try:
            names = survey_names(surveys, ("gaia", "jwst", "roman", "roman-sim"))
            return cone(ra, dec, radius_arcsec=radius_arcsec, limit=limit,
                        surveys=names, program=roman_program)
        except ValueError as error:
            raise HTTPException(422, detail=str(error)) from error

    @app.get("/v1/gaia/{source_id}/photometry")
    def gaia_photometry(source_id: str):
        return archive_call("gaia", gaia_epoch_photometry, source_id)

    @app.get("/v1/jwst/{obsid}/products")
    def jwst_product_list(obsid: str, download: bool = False, directory: str = "data/jwst"):
        return archive_call("jwst", jwst_products, obsid, download=download, directory=directory)

    @app.post("/v1/experiments/orbital-decay")
    def run_experiment(request: RunRequest):
        if not partitions:
            raise HTTPException(409, detail="Start sky serve with --ztf to enable runs.")
        try:
            context = RunContext(partitions, survey_names(request.surveys),
                                 top=request.top, workers=workers,
                                 limit_per_partition=request.limit_per_partition)
            return orbital_decay(context)
        except ValueError as error:
            raise HTTPException(422, detail=str(error)) from error

    return app
