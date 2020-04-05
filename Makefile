
DOCKERFILES = $(wildcard *.Dockerfile)
BASE_REPO = hpobiolab
date = $(shell date +"%Y-%b-%d")

SINGULARITY_DIR=singularity

make_tags:
	mkdir -p make_tags

push: make_tags/$(APP).built make_tags


build-cacheless: dockerfiles/$(APP).Dockerfile make_tags
	+docker build --no-cache -t $(BASE_REPO)/$(APP) -f dockerfiles/$(APP).Dockerfile . && touch make_tags/$(APP).built
	+docker tag $(BASE_REPO)/$(APP) $(BASE_REPO)/$(APP):latest
	+docker tag $(BASE_REPO)/$(APP) $(BASE_REPO)/$(APP):$(date)


push-cacheless: build-cacheless  make_tags
	+docker push $(BASE_REPO)/$(APP):latest && touch make_tags/$(APP).latest
	+docker push $(BASE_REPO)/$(APP):$(date) && touch make_tags/$(APP).$(date)

singularity: push
	+docker run quay.io/singularity/docker2singularity $(ARGS) --name $(BASE_REPO)/$(APP) $(BASE_REPO)/$(APP)



build: dockerfiles/$(APP).Dockerfile make_tags
	+docker build -t $(BASE_REPO)/$(APP) -f dockerfiles/$(APP).Dockerfile . && touch make_tags/$(APP).built
	+docker tag $(BASE_REPO)/$(APP) $(BASE_REPO)/$(APP):latest
	+docker tag $(BASE_REPO)/$(APP) $(BASE_REPO)/$(APP):$(date)

push: build  make_tags
	+docker push $(BASE_REPO)/$(APP):latest && touch make_tags/$(APP).latest
	+docker push $(BASE_REPO)/$(APP):$(date) && touch make_tags/$(APP).$(date)

base: base/base.Dockerfile
	docker build -t erictdawson/base -f base/base.Dockerfile . && touch make_tags/base.built
	+docker tag erictdawson/base erictdawson/base:latest
	+docker tag erictdawson/base erictdawson/base:$(date)

clean:
	$(RM) make_tags

.PHONY: build push clean
