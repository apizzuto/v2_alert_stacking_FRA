import setuptools

long_message = 'FRANCIS: Fast Response Analysis for Neutrino Coincidences with IceCube Signals'
version = "0.0.1"

setuptools.setup(
    name="francis", 
    version=version,
    author="Pizzuto, Alex",
    author_email="",
    description="Code for following up icecube alert events",
    long_description=long_message,
    #long_description_content_type="text/markdown",
    url="https://github.com/apizzuto/v2_alert_stacking_FRA",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
    ],
    python_requires='>=3.1',
)
