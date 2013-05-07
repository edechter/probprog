

(declare (usual-integrations))

(load "pp")
(load "pp-interface")
(load "pp-records")

(load "generic-rv")
(load "beta")
(load "bern")

;; a global beta-bernoulli data structure that keeps track of
;; bernoullis and betas
;; KEYS: beta random variables
;; VALUES: objects of record type beta-bern-ds

(define *betas*)
(define (*betas*:init!)
  (set! *betas* (make-eq-hash-table)))
(*betas*:init!)
  
(define-record-type <beta-bern-ds>
  (beta-bern-ds:new a b berns)
  beta-bern-ds?
  (a beta-bern-ds:a beta-bern-ds:set-a!)
  (b beta-bern-ds:b beta-bern-ds:set-b!)
  (berns beta-bern-ds:berns beta-bern-ds:set-berns!))

(define (*betas*:init-new! beta-rv a b)
  (let ((ds (beta-bern-ds:new a b (make-eq-hash-table))))
    (hash-table/put! *betas* beta-rv ds)
    ds))

(define (*betas*:remove-beta! beta-rv)
  (hash-table/remove! *betas* beta-rv))

(define (*betas*:get-beta beta-rv)
  (hash-table/get *betas* beta-rv 'no-value))

(define (*betas*:add-bern! beta-rv bern-rv val)
  (let* ((berns (beta-bern-ds:berns (*betas*:get-beta beta-rv))))
         (hash-table/put! berns bern-rv val)))

(define (beta-bern-ds:show ds)
  (display "a: ") (display (beta-bern-ds:a ds)) (newline)
  (display "b: ") (display (beta-bern-ds:b ds)) (newline)
  (display "berns: ") (display (hash-table->alist (beta-bern-ds:berns ds))) (newline))

(define (*betas*:show)
  (for-each beta-bern-ds:show
            (hash-table/datum-list *betas*)))

(define (beta-bern-ds:counts ds)
  (let* ( (counts (hash-table/datum-list (beta-bern-ds:berns ds)))
          (heads (filter (lambda (x) (eq? 1 x)) counts))
          (tails (filter (lambda (x) (eq? 0 x)) counts))
          (new-a (+ (length heads) (beta-bern-ds:a ds)))
          (new-b (+ (length heads) (beta-bern-ds:b ds))))
    (cons a b)))


(define (*betas*:counts beta-rv)
  (beta-bern-ds:counts (*betas*:get-beta beta-rv)))          

(define (beta a b #!optional proposer)
  (if (default-object? proposer)
      (set! proposer (proposals:logistic-gaussian 
                      (/ (beta:var (beta:make-params  a b))  4))))
  
  (let ( (rv (rv:new 'beta
                     the-unsampled-value
                     (pair->flo-vector (cons a b))
                     beta:rvs
                     beta:log-likelihood
                     beta:force-hook
                     beta:force-set-hook
                     proposer)))
    (*betas*:init-new! rv a b)
    rv))

(define ((beta:rvs-global beta-rv) params)
  (let ((a&b (*betas*:counts beta-rv)))
    (beta-rvs (pair->flo-vector a&b))))

(define (beta:force-hook rv)
  'unspecific)

(define (beta:force-set-hook rv)
  'unspecific)

(define (bern p #!optional proposer)
  (if (default-object? proposer)
      (set! proposer (proposals:flip 0.5)))
  (cond ( (number? p)
          (let ((rv (rv:new 'bern
                            the-unsampled-value
                            p
                            bern:rvs 
                            bern:log-likelihood 
                            'unspecific
                            'unspecific
                            proposer)))
            rv))
          ;; if p is an unsampled beta rv
        ((and (rv:type? p 'beta) (rv:unsampled? p))
         (let ((rv (rv:new 'bern
                           the-unsampled-value
                           ;; params
                           p
                           ;; sampler
                           (lambda (beta-rv)
                               (bern:rvs (rv:force beta-rv)))
                           ;; marginal log-likelihood TODO: currently
                           ;; assumes there is only one beta in the
                           ;; world
                           (lambda (x beta-rv)
                             (let ((ds (*betas*:get beta-rv)))
                               (bern:marginal-log-likelihood x
                                                             (beta-bern-ds:a ds)
                                                             (beta-bern-ds:b ds))))

                           ;; force-hook
                           (bern:force-hook p)
                           (bern:force-set-hook p)
                           proposer)))
             rv))
          ((and (rv:type? p 'beta) (rv:sampled? p))
           (let ((rv (rv:new 'bern
                           the-unsampled-value
                           ;; params
                           p
                           ;; sampler
                           (lambda (beta-rv)
                               (bern:rvs (rv:val beta-rv)))
                           ;; conditional log-likelihood 
                           (lambda (x beta-rv)
                             (bern:rvs x (rv:val beta-rv)))
                           ;; force-hook
                           (bern:force-hook p)
                           ;; force-set-hook
                           (bern:force-set-hook p)
                           proposer)))
             rv))
          ( (error "BERN: don't know what to do with this type of parameter." p))))

(define ((bern:force-hook beta-rv) bern-rv)
  (*betas*:modify-bern! beta-rv bern-rv (rv:val bern-rv)))
  
  
(define ((bern:force-set-hook beta-rv) bern-rv )
  (display "Adding bernoulli to beta!") (newline)
  (*betas*:add-bern! beta-rv bern-rv (rv:val bern-rv)))

(define (emit-bern bern val)
  (force-set! bern val))
  
;; Goal program
(define (test-beta-bern)
  (let* ((p (beta 1 2))
        (x1 (bern p))
        (x2 (bern p))
        (x3 (bern p))
        (x4 (bern p))
        (x5 (bern p))
        (x6 (bern p)))
    (emit-bern x1 0)
    (emit-bern x2 0)
    (emit-bern x3 0)
    (emit-bern x4 0)
    (emit-bern x5 0)
    (emit-bern x6 0)
    p))





  


