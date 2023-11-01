FROM debian:bookworm-slim AS build
COPY . /app
WORKDIR /app
RUN apt-get update && apt-get install -y build-essential && apt-get clean
RUN make release

FROM debian:bookworm-slim
COPY --from=build /app/bin/dpsim /app/
RUN apt-get update && apt-get install -y gnuplot-nox imagemagick && apt-get clean
CMD ["/app/dpsim"] 
