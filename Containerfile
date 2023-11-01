FROM debian:bookworm-slim AS build
RUN apt-get update && apt-get install -y build-essential && apt-get clean
COPY . /app
WORKDIR /app
RUN make release

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y gnuplot-nox imagemagick && apt-get clean
COPY --from=build /app/bin/dpsim /app/
CMD ["/app/dpsim"] 
